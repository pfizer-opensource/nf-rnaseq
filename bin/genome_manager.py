#!/usr/bin/env python3

import os
import fileinput
from pathlib import Path
import argparse
import logging
import grp
import json
from platform import python_version
import shutil
import yaml
from pydantic import BaseModel, validator, FilePath, DirectoryPath, ValidationError
from typing import Optional, Union
import hashlib
from getpass import getuser
from generate_gtf_entry import YamlGeneCollection, YamlGeneModel
import generate_gtf_entry

__author__ = 'Rob Moccia'
__version__ = '0.1'


class DuplicateGenomeError(Exception):
    pass
class FileMismatchError(Exception):
    pass

class FileFormatError(Exception):
    pass

class GenomeFile(BaseModel):
    """Schema for a file that is part of a genome."""

    type: str
    path: dict[str, FilePath] # one each for e.g., hpc, aws, gpfs
    checksum: str=None
    source: Optional[str]
    content: Optional[str]
    parent: Optional[FilePath]
    
    @validator('type', pre=True)
    def validate_type(cls, value):
        valid_types = ['fasta', 'gtf', 'refflat', 'rrna_interval_list', 'gff3', 'yaml_gene_model']
        if value not in valid_types:
            raise ValueError(f'{value} is not a recognized type ({valid_types})')
        return value

    @validator('checksum', always=True)
    def add_checksum(cls, val, values):
        target = next(iter(values['path'].values()))
        filesize = Path(target).stat().st_size
        if filesize > 100000:
            return filesize
        with open(target, 'rb') as f:
            file_hash = hashlib.md5()
            while chunk := f.read(16384):
                file_hash.update(chunk)
        return file_hash.hexdigest()

    @validator('source')
    def validate_source(cls, value):
        valid_sources = ['genome', 'transcriptome']
        if value is None:
            return
        if value not in valid_sources:
            raise ValueError(f'{value} is not a recognized source ({valid_sources})')
        return value

    # if the file is specified as type yaml_gene_model must validate that it
    # follows the specification for YamlGeneModel or YamlGeneCollection in generate_gtf_entry.py
    @validator('path', pre=True)
    def validate_path(cls, value, values):
        if values['type'] == 'yaml_gene_model':
            for system_name, yaml_file in value.items():
                with open(yaml_file, 'r') as f:
                    try:
                        input_obj = yaml.load(f, Loader=yaml.CLoader)
                    except:
                        logger.exception(
                            f'type = yaml_gene_model and {yaml_file} (system_name: {system_name}) is not a valid YAML file')
                        raise
                    try:
                        if isinstance(input_obj, dict):
                            YamlGeneCollection.parse_obj(input_obj)
                        elif isinstance(input_obj, list):
                            YamlGeneModel.parse_obj(input_obj[0])
                        else:
                            raise ValidationError
                    except ValidationError:
                        logger.exception(
                            f'type = yaml_gene_model and {yaml_file} ({system_name}) is not a valid YAML gene model specification')
                        raise
            return value
        else:
            return value
class GenomePath(BaseModel):
    """Schema for a directory that is part of a genome."""

    path: dict[str, DirectoryPath] # one each for e.g., hpc, aws, gpfs
    type: str
    source: Optional[str]
    content: Optional[str]
    parent: Optional[FilePath]
    
    @validator('type')
    def validate_type(cls, value):
        valid_types = ['star_index']
        if value not in valid_types:
            raise ValueError(f'{value} is not a recognized type ({valid_types})')
        return value

    @validator('source')
    def validate_source(cls, value):
        valid_sources = ['genome', 'transcriptome']
        if value is None:
            return
        if value not in valid_sources:
            raise ValueError(f'{value} is not a recognized source ({valid_sources})')
        return value

class BaseGenome(BaseModel):
    """Schema for the core files representing a genome assembly (i.e., genome fasta and GTF)"""

    id: str
    species: str
    version: int
    assembly: str
    assembly_type: str
    sequence_type: str
    genome_fasta: GenomeFile
    gtf: GenomeFile
    description: Optional[str]

    # @validator('assembly_type')
    # def validate_assembly_type(cls, value):
    #     valid_assembly_types = ['pa', 'tl', 'primary_assembly', 'toplevel']
    #     if value not in valid_assembly_types:
    #         raise ValueError(f'{value} is not a valid assembly type ({valid_assembly_types})')
    #     return value

    @validator('sequence_type')
    def validate_sequence_type(cls, value):
        valid_sequence_types = ['dna', 'dna_rm', 'dna_sm']
        if value not in valid_sequence_types:
            raise ValueError(f'{value} is not a valid sequence type ({valid_sequence_types})')
        return value

class Genome(BaseModel):
    """Schema for representing a genome assembly, annotation, and associated files and metadata"""

    id: str
    base: BaseGenome
    transcriptome_fasta: GenomeFile
    star_index: GenomePath
    refflat: GenomeFile
    rrna_interval_list: GenomeFile
    description: Optional[str]

class GenomeCollection(BaseModel):
    """Schema for representing a collection of Genome objects"""

    genomes: dict[str, Genome]

class UserDefinedGene(BaseModel):
    """Schema for representing a transcript(s) that can be added to a genome"""
    
    gene_model: dict[int, GenomeFile]
    fasta: GenomeFile
    id: str

    class Config:
        validate_assignment = True

    @validator('gene_model')
    def validate_gene_model(cls, val: dict[int, GenomeFile]):
        """Enforce the same gene_id in every YAML associated with the same UserDefinedGene"""
        ids = set()
        for genome_file in val.values():
            for _,yaml_file in genome_file.path.items():
                yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
                current_id = YamlGeneModel.parse_obj(yaml_obj[0]).gene_id
                ids.add(current_id)
                if len(ids) > 1:
                    raise ValueError(f'YAML files for different versions have different gene_id: {ids}')
        return val

    @validator('fasta')
    def validate_fasta(cls, val: GenomeFile, values) -> GenomeFile:
        """Check that fasta file has only one entry named exactly the same as gene_id in the YAML model"""
        check_file = next(iter(values['gene_model'].values()))
        yaml_file = next(iter(check_file.path.values()))
        yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
        gene_id = YamlGeneModel.parse_obj(yaml_obj[0]).gene_id
        for system_name, filename in val.path.items():
            description = []
            with open(filename, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        continue
                    else:
                        description.append(line)
            if len(description) != 1:
                logger.error(
                    f'fasta file {filename} for system_name {system_name} appears to have more than one sequence entry')
                raise ValueError('fasta file must have only 1 entry')
            else:
                fasta_id = description.pop().strip().lstrip('>')
                if fasta_id != gene_id:
                    logger.error(
                        f'sequence name in {filename} ({fasta_id}) does not match gene_id in YAML ({gene_id})')
                    raise ValueError('fasta file sequence name must match gene_id in YAML gene model')
        return val

    @validator('id')
    def validate_id(cls, val: str, values) -> str:
        check_file = next(iter(values['gene_model'].values()))
        yaml_file = next(iter(check_file.path.values()))
        yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
        if val != (gene_id := YamlGeneModel.parse_obj(yaml_obj[0]).gene_id):
            raise ValueError(f'provided id ({val}) does not match gene_id in {yaml_file} ({gene_id})')
        return val

    def get_version(self, version: int, system_name: str) -> Union[str, bytes, os.PathLike]:
        """
        Return path to a gene_model version by version number and system name.
        Latest version is represented by code -1.
        """
        version = int(version)
        if version < 0:
            version = sorted(self.gene_model.keys())[-1]
        return Path(self.gene_model[version].path[system_name])

    def add_version(self, yaml_file: Union[str, bytes, os.PathLike],
            system_name: str) -> tuple[int, Union[str, bytes, os.PathLike]]:
        """
        Add a new YAML gene_model to the gene_model dict, auto-incrementing the version
        Returns a tuple of the newly assigned version number as well as path where YAML
        was written which can be used to restore registry to pre-existing state on downstream
        failure.
        """
        latest_version = sorted(self.gene_model.keys())[-1]
        this_version = latest_version + 1
        # ensure consistent path for system_name by recovering it from a previous version
        try:
            prior_yaml_path = next(iter(self.gene_model.values())).path[system_name]
            registry_path = prior_yaml_path.parent
        except KeyError:
            logger.exception(f'{system_name} is not a valid system_name for {self.id}')
            raise
        except:
            logger.exception(f'failed to add new version to {self.id}')
            raise
        yaml_dest = Path(registry_path,
            Path(yaml_file).stem + '_v' + str(this_version).zfill(2) + Path(yaml_file).suffix)
        try:
            shutil.copy(yaml_file, yaml_dest)
            new_model = GenomeFile(path={system_name: yaml_dest}, type='yaml_gene_model')
            for version, next_model in self.gene_model.items():
                if new_model.checksum == next_model.checksum:
                    logger.error(f'checksum matches model stored in version {version}')
                    raise ValueError(f'YAML gene model is identical to a previously stored version')
            self.gene_model[latest_version + 1] = new_model
        except:
            logger.error(f'failed to add {yaml_file} as new version for {self.id}')
            if yaml_dest.exists():
                yaml_dest.unlink()
                logger.info(f'ERROR RECOVERY: undid addition of {yaml_dest} to registry')
            else:
                logger.info(
                    'ERROR RECOVERY: no files were written so no action necessary to restore registry to pre-existing state')
            raise
        return latest_version + 1, yaml_dest

class UserDefinedGeneCollection(BaseModel):
    """Schema for representing all available modifications that can be added to a genome"""

    modifications: dict[str, UserDefinedGene]


def initialize(registry_path: Union[str, bytes, os.PathLike], group_name: Optional[str]=None, **kwargs) -> None:
    """Initialize a new genome registry"""
    os.makedirs(registry_path)

    open_subdirs = ['user_defined_genes', '.user-registry', '.log']
    restricted_subdirs = ['genomes', '.genome-registry']
    for dirname in open_subdirs:
        os.makedirs(os.path.join(registry_path, dirname))
    for dirname in restricted_subdirs:
        os.makedirs(os.path.join(registry_path, dirname))
        os.chmod(os.path.join(registry_path, dirname), mode=0o775)

    if group_name:
        uid = os.getuid()
        gid = grp.getgrnam(group_name).gr_gid
        os.chown(registry_path, uid, gid)
        for dirname in open_subdirs + restricted_subdirs:
            os.chown(os.path.join(registry_path, dirname), uid, gid)

def gene_model_from_yaml(yaml_file: Union[str, bytes, os.PathLike]) -> YamlGeneModel:
    """Create YamlGeneModel object from a YAML file"""
    try:
        yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
    except:
        raise
    if isinstance(yaml_obj, dict):
        raise ValueError(f'{yaml_file} has top-level key(s) -- consider register-gene-collection')
    elif not isinstance(yaml_obj, list):
        raise ValueError(f'{yaml_file} is not a valid YAML gene model')
    else:
        if len(yaml_obj) > 1:
            raise ValueError(f'{yaml_file} contains {len(yaml_obj)} entries -- remove extra entries or re-format as collection')
        model = YamlGeneModel.parse_obj(yaml_obj[0])
        return model

def validate_user_gene_file(filename:  Union[str, bytes, os.PathLike]) -> None:
    """
    Validate that a file ends with a newline. This is necessary for files that will be
    concatenated prior to returning like those for user-defined genes.
    """
    with open(filename, 'r') as f:
        filestring = f.read()
    if not filestring.endswith('\n'):
        logger.error(f'{filename} must end with a newline character')
        raise FileFormatError(f'{filename} does not end with a newline')

def build_new_user_defined_gene(fasta: Union[str, bytes, os.PathLike], yaml_file: Union[str, bytes, os.PathLike],
        system_name: str, registry_path: Union[str, bytes, os.PathLike], **kwargs) -> UserDefinedGene:
    """
    Create a UserDefinedGene object from fasta file and YAML gene model for a single gene.
    If trying to add a collection of genes (e.g., ERCC) as a multi-fasta and associated YAML
    gene model definitions see build_new_user_defined_gene_collection().
    """
    try:
        model = gene_model_from_yaml(yaml_file)
    except:
        logger.exception('register-gene failed while loading gene model from YAML')
        raise
    # always use absolute path in registry entries
    registry_path = Path(registry_path).resolve()
    target_dir = Path(registry_path, 'user_defined_genes', model.gene_id)
    if target_dir.exists():
        logger.error(f'{model.gene_id} has already been added to registry ({target_dir})')
        raise FileExistsError(f'{target_dir} already exists  -- try update-yaml-model')
    # fasta_dest = Path(target_dir, fasta)
    fasta_dest = Path(target_dir, f'{model.gene_id}.fa')
    # yaml_fname = Path(Path(yaml_file).stem + '_v01' + Path(yaml_file).suffix)
    yaml_fname = f'{model.gene_id}_v01.yaml'
    yaml_dest = Path(target_dir, yaml_fname)
    try:
        target_dir.mkdir()
        shutil.copy(fasta, fasta_dest)
        shutil.copy(yaml_file, yaml_dest)
        yaml_file = GenomeFile(path={system_name: yaml_dest}, type='yaml_gene_model')
        fasta_file = GenomeFile(path={system_name: fasta_dest}, type='fasta')
        new_gene = UserDefinedGene(id=model.gene_id, fasta=fasta_file, gene_model={1: yaml_file})
        logger.info(f'built new user-defined-gene {model.gene_id}')
    except:
        logger.exception(f'failed to build gene for fasta: {fasta}, YAML: {yaml_file}')
        logger.info('starting error recovery')
        try:
            shutil.rmtree(target_dir)
            logger.info(f'ERROR RECOVERY: {target_dir} successfully removed')
        except FileNotFoundError:
            logger.info(f'ERROR RECOVERY: {target_dir} was not created - no action required')
            raise
        raise
    return new_gene

def register_user_defined_gene(fasta: Union[str, bytes, os.PathLike], yaml_file: Union[str, bytes, os.PathLike],
        system_name: str, registry_path: Union[str, bytes, os.PathLike], **kwargs) -> None:
    """
    Add a new user-defined gene to the registry. Called by command line via `register-gene`
    See also update_user_defined_gene() to add a new gene model to an existing user-defined gene.
    """
    logger.info(f'{getuser()} called register-gene for registry: {Path(registry_path).resolve()}')
    validate_user_gene_file(fasta)
    validate_user_gene_file(yaml_file)

    new_gene = build_new_user_defined_gene(fasta=fasta, yaml_file=yaml_file, system_name=system_name,
        registry_path=registry_path)
    registry_file = Path(registry_path, '.user-registry', new_gene.id + '.conf')
    with registry_file.open('w') as f:
        f.write(new_gene.json())
        logger.info(f'{getuser()} added user-defined gene {new_gene.id} version 1 to registry')

def update_user_defined_gene(registry_path: Union[str, bytes, os.PathLike],
        yaml_file: Union[str, bytes, os.PathLike], system_name: str, **kwargs) -> None:
    """
    Adds a new gene model to an existing user-defined gene that was added via `register-gene`
    Function is called by command line argument `update-gene`
    """
    logger.info(f'{getuser()} called update-gene for system_name: {system_name}; yaml_file: {yaml_file}')
    validate_user_gene_file(yaml_file)

    # load yaml_file and parse with YamlGeneModel
    try:
        model = gene_model_from_yaml(yaml_file)
    except:
        logger.exception('update-gene failed while loading gene model from YAML')
        raise
    # find the .conf file
    # if it doesn't exist, then this isn't an update operation so fail with error message
    registry_path = Path(registry_path).resolve()
    registry_file = Path(registry_path, '.user-registry', model.gene_id + '.conf')
    if not registry_file.exists():
        logger.error(
            f'update-gene failed: no configuration file for {model.gene_id} found at {registry_file}')
        raise FileNotFoundError(f'{registry_file.resolve()} not found; try register-gene if adding gene for the first time')
    gene = load_user_defined_gene(registry_file)

    # copy YAML to registry and update JSON config
    new_version_num, yaml_dest = gene.add_version(yaml_file, system_name)
    # UserDefinedGene.add_version() has error handling implemented so if it returns,
    # addition was successful and it's safe to update the config file.
    # However, an error writing the new config file would put the registry in a corrupted state.
    # To be safe, hold a copy of the original that can be restored and delete the added YAML file
    # on failure at this step.
    try:
        # hold onto a backup copy of original config in case it needs to be restored on failure
        with registry_file.open('r') as f:
            original_registry_file = f.read()
        with registry_file.open('w') as f:
            f.write(gene.json())
            logger.info(f'{getuser()} successfully updated {gene.id} to version {new_version_num}')
    except:
        logger.error(f'update-gene for {model.gene_id} encountered error updating the registry config file')
        logger.info(f'ERROR RECOVERY: restoring previous config file version')
        logger.info(f'{model.gene_id} config JSON provided here as failsafe')
        logger.info(f'{original_registry_file}')
        with registry_file.open('w') as f:
            f.write(original_registry_file)
        if Path(yaml_dest).exists():
            logger.info(f'ERROR RECOVERY: deleting {yaml_dest} from registry')
            Path(yaml_dest).unlink()
        else:
            logger.info(f'ERROR RECOVERY: {yaml_file} was not written to registy - no further action required')

def load_user_defined_gene(registry_file: Union[str, bytes, os.PathLike]) -> UserDefinedGene:
    """Load a UserDefinedGene object from the registry"""
    try:
        with open(registry_file, 'r') as f:
            json_obj = json.load(f)
        logger.info(f'loaded {registry_file}')
        gene = UserDefinedGene.parse_obj(json_obj)
        logger.info(f'parsed {gene.id} gene model from {Path(registry_file).resolve()}')
    except:
        logger.exception(f'failed to load gene model from {Path(registry_file).resolve()}')
        raise
    return gene

def get_user_defined_genes(registry_path: Union[str, bytes, os.PathLike], gene_ids: Union[str, list[str]],
        system_name: str, outdir: Union[str, bytes, os.PathLike] = Path('.'), 
        version_delim: str='.', **kwargs) -> None:
    """
    Retrieve user-defined gene(s) and write files to disk. If multiple genes, concatenate files
    before returning.
    Called from command line with command `get-genes`
    """
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]
    if not isinstance(gene_ids, list):
        raise TypeError('gene_ids must be of type str or list of strings')
    if not all([isinstance(gene, str) for gene in gene_ids]):
        raise TypeError('gene_ids must be a string or list of strings')

    fasta_files = []
    yaml_files = []
    collected_ids = []

    for gid in gene_ids:
        if version_delim in gid:
            gene_id, version = gid.strip().split(version_delim)
        else:
            gene_id = gid.strip()
            version = -1

        registry_file = Path(registry_path, '.user-registry', gene_id + '.conf')
        gene = load_user_defined_gene(registry_file)
        fasta_files.append(Path(gene.fasta.path[system_name]))
        yaml_files.append(gene.get_version(version, system_name))
        collected_ids.append(gene.id)
        print(collected_ids)
        print(gene_id)
        print(version)
    if len(collected_ids) > 3:
        basename = 'custom'
    else:
        basename = '.'.join(collected_ids)

    with open(Path(outdir, basename + '.fa'), 'wb') as outfile:
        for f in fasta_files:
            with open(f, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

    yaml_concat = fileinput.input(yaml_files)
    yaml_out = Path(outdir, basename + '.yaml')
    with open(yaml_out, 'w') as outfile:
        for line in yaml_concat:
            outfile.write(line)

    gene_model = generate_gtf_entry.read_genes_from_yaml(yaml_out)
    gtf = generate_gtf_entry.generate_gtf(gene_model)
    with open(Path(outdir, basename + '.gtf'), 'w') as outfile:
        outfile.write(gtf)

def load_genome(registry_file: Union[str, bytes, os.PathLike]) -> GenomeCollection:
    """Parse a GenomeCollection object from a registry JSON filepath"""
    try:
        model = GenomeCollection.parse_file(Path(registry_file))
        logger.info(f'loaded GenomeCollection from {Path(registry_file).resolve()}')
        return model
    except:
        logger.exception(f'failed to load from {Path(registry_file).resolve()}')
        raise

def copy_with_logging(src: Union[str, bytes, os.PathLike], dest: Union[str, bytes, os.PathLike]
        ) ->  Union[str, bytes, os.PathLike]:
    """Wrapper around shutil.copy and shutil.copytree that adds logging"""
    src = Path(src)
    dest = Path(dest).resolve()
    if src.is_file():
        try:
            logger.info(f'copying {src} to {dest}')
            if not dest.parent.exists():
                dest.parent.mkdir(parents=True)
            if not dest.exists():
                write_path = shutil.copy(src, dest)
                logger.info(f'successfully copied {src} to {dest}')
            else:
                logger.error(f'attempted to overwrite existing file {dest} when copying {src}')
                raise FileExistsError(f'{dest} already exists')
        except:
            logger.exception(f'failed to copy {src} to {dest}')
            raise
    elif src.is_dir():
        try:
            logger.info(f'copying directory tree at {src} to {dest}')
            write_path = shutil.copytree(src, dest)
            logger.info(f'successfully copied directory tree at {src} to {dest}')
        except:
            logger.exception(f'failed to copy directory tree at {src} to {dest}')
            raise
    else:
        logger.exception(f'failed to copy {src} to {dest} because source is not a file or a directory')
        raise TypeError(f'{src} is not a file or a directory')
    return write_path

def write_genome_files(genome_dir: Union[str, bytes, os.PathLike], genome_fasta: Union[str, bytes, os.PathLike],
        gtf: Union[str, bytes, os.PathLike], transcriptome_fasta: Union[str, bytes, os.PathLike],
        refflat: Union[str, bytes, os.PathLike], rrna_interval_list: Union[str, bytes, os.PathLike],
        star_index: Union[str, bytes, os.PathLike]) -> dict[str, Union[str, bytes, os.PathLike]]:
    """
    Write genome files to the registry enforcing a consistent directory structure.
    Currently implemented specifically for Ensembl genomes.
    """
    genome_fasta_path = copy_with_logging(genome_fasta, Path(genome_dir, 'source', genome_fasta.name))
    gtf_path = copy_with_logging(gtf, Path(genome_dir, 'source', gtf.name))
    transcriptome_fasta_path = copy_with_logging(transcriptome_fasta, Path(genome_dir, 'derived', transcriptome_fasta.name))
    refflat_path = copy_with_logging(refflat, Path(genome_dir, 'derived', refflat.name))
    rrna_interval_list_path = copy_with_logging(rrna_interval_list, Path(genome_dir, 'derived', rrna_interval_list.name))
    star_index_path = copy_with_logging(star_index, Path(genome_dir, star_index.name))
    return {
        'genome_fasta': genome_fasta_path,
        'gtf': gtf_path,
        'transcriptome_fasta': transcriptome_fasta_path,
        'refflat': refflat_path,
        'rrna_interval_list': rrna_interval_list_path,
        'star_index': star_index_path
        }

def build_new_genome(params: dict, registry_path: Union[str, bytes, os.PathLike]) -> Genome:
    """
    Create a new Genome object using a dictionary of parameters.
    Format requires top level keys `base` (BaseGenome params), `genome_fasta`,
    `gtf`, `transcriptome_fasta`, `star_index`, `refflat`, and `rrna`.
    """
    # build registry directory name for this release and species
    genome_dir = Path(registry_path, 'genomes', f"release-{params['base']['version']}",
        f"{params['base']['species']}")
    if genome_dir.exists():
        logger.error(f'a genome has already been stored in {genome_dir.resolve()}')
        raise FileExistsError(f'a genome has already been stored in {genome_dir.resolve()}')
    try:
        write_paths = write_genome_files(
            genome_dir=genome_dir,
            genome_fasta=Path(params['genome_fasta']),
            gtf=Path(params['gtf']),
            transcriptome_fasta=Path(params['transcriptome_fasta']),
            refflat=Path(params['refflat']),
            rrna_interval_list=Path(params['rrna_interval_list']),
            star_index=Path(params['star_index'])
        )
    except FileExistsError:
        logger.exception('aborting due to attempted overwrite of existing registry files')
        raise

    try:
        # the paths returned in write_paths are already absolute paths
        system_name = params['system_name']
        genome_fasta = GenomeFile(type='fasta', source='genome',
            path={system_name: write_paths['genome_fasta']})
        gtf = GenomeFile(type='gtf', path={system_name: write_paths['gtf']})
        transcriptome_fasta = GenomeFile(type='fasta', source='transcriptome',
            path={system_name: write_paths['transcriptome_fasta']})
        star_index = GenomePath(type='star_index', path={system_name: write_paths['star_index']})
        refflat = GenomeFile(type='refflat', path={system_name: write_paths['refflat']})
        rrna_interval_list = GenomeFile(type='rrna_interval_list', path={system_name: write_paths['rrna_interval_list']})
        base_genome = BaseGenome(
            gtf=gtf,
            genome_fasta=genome_fasta,
            **params['base'])
        new_genome = Genome(
            id=base_genome.id,
            base=base_genome,
            transcriptome_fasta=transcriptome_fasta,
            star_index=star_index,
            refflat=refflat,
            rrna_interval_list=rrna_interval_list)
        logger.info(f'successfully built new genome {new_genome.id}')
        return new_genome
    except:
        logger.exception('Exception raised in build_new_genome(); starting ERROR RECOVERY')
        # try to delete any files and directories that were successfully created
        try:
            shutil.rmtree(genome_dir)
            logger.info(f'ERROR RECOVERY: {genome_dir} successfully removed')
        except FileNotFoundError:
            logger.info(f'ERROR RECOVERY: {genome_dir} was not created - no action required')
        # if this is the only species in this release, remove the release directory too
        if genome_dir.parent.is_dir():
            try:
                genome_dir.parent.rmdir()
                logger.info(f'ERROR RECOVERY: removed empty release directory {genome_dir.parent}')
            except:
                logger.exception(f'failed to remove {genome_dir.parent}')
                raise
        raise

def parse_json_genome_args(json_args: Union[str, bytes, os.PathLike]) -> dict:
    """Read genome args from JSON file and return as dictionary with id converted to lowercase"""
    try:
        with open(json_args, 'r') as f:
            params = json.load(f)
        # convert id to lower by convention to enable case-insensitive searching
        params['base']['id'] = params['base']['id'].lower()
    except FileNotFoundError:
        logger.exception(f"failed to load {Path(json_args).resolve()}")
        raise
    except KeyError:
        logger.exception(f'{Path(json_args).resolve()} is not formatted properly')
        raise
    except AttributeError:
        logger.exception('failed to convert id to lowercase')
        raise
    except:
        logger.exception()
        raise
    return params

def register_genome(json_args: str, registry_path: str, system_name: str, **kwargs) -> None:
    """
    Add a new Genome to a GenomeCollection and write the json configuration file.
    Called via command line by `register-genome`.
    """
    registry_path = Path(registry_path).resolve()
    params = parse_json_genome_args(json_args)
    params['system_name'] = system_name
    logger.info(f"{getuser()} called register-genome for {params['base']['id']}")
    conf_filename = f"{params['base']['version']}.json"
    registry_file = Path(registry_path, '.genome-registry', conf_filename)
    if registry_file.exists():
        genome_registry = load_genome(registry_file)
        if params['base']['id'] in genome_registry.genomes.keys():
            logger.error(f"aborting - genome with id {params['base']['id']} was already registered")
            raise DuplicateGenomeError(f"there is already a genome with id: {params['base']['id']}")
        else:
            genomes_list = genome_registry.genomes
    else:
        genomes_list = dict()
        # registry_file.parent.mkdir(exist_ok=True)

    try:
        genome = build_new_genome(params, registry_path)
        genomes_list[genome.id] = genome
        with registry_file.open('w') as f:
            f.write(GenomeCollection(genomes=genomes_list).json())
            logger.info(f'{getuser()} added genome {genome.id} to registry')
    except:
        logger.exception('exception raised in register_genome()')
        raise

def delete_genome(genome):
    """Remove a genome from the registry"""
    pass

def start_logger(registry_path: Union[str, bytes, os.PathLike], **kwargs) -> None:
    if kwargs['command'] == 'init':
        # no logging for init
        return
    elif kwargs['command'] in ['register-genome']:
        logfilename = 'genome-manager.log'
    elif kwargs['command'] in ['register-gene', 'update-gene']:
        basename = os.path.splitext(os.path.basename(kwargs['yaml_file']))[0]
        logfilename = basename + '.log'
    elif kwargs['command'] == 'get-genes':
        logfilename = 'get-genes.log'
    else:
        raise NotImplementedError(f"logger not implemented for {kwargs['command']}")
    logging.basicConfig(
        filename=os.path.join(registry_path, '.log', logfilename),
        encoding='utf-8',
        level=logging.INFO,
        format='%(asctime)s %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f'genome_manager.py v{__version__} called by {getuser()} using Python {python_version()}')
    return logging.getLogger()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))
    sp = parser.add_subparsers(dest='command')

    register_genome_parser = sp.add_parser('register-genome', help='register a genome')
    register_genome_parser.set_defaults(func=register_genome)
    register_genome_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    register_genome_parser.add_argument('--json-args', required=True,
        help = 'path to JSON file containing function arguments')
    register_genome_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., hpc-durham, AWS, etc.)')

    register_gene_parser = sp.add_parser('register-gene', help='register a user-defined gene')
    register_gene_parser.set_defaults(func=register_user_defined_gene)
    register_gene_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    register_gene_parser.add_argument('--yaml-file', required=True,
        help = 'path to YAML file containing gene model for custom fasta')
    register_gene_parser.add_argument('--fasta', required=True,
        help = 'path to fasta file containing sequence of a genome modification')
    register_gene_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., hpc-durham, AWS, etc.)')

    init_parser = sp.add_parser('init', help='initialize a new genome registry')
    init_parser.set_defaults(func=initialize)
    init_parser.add_argument('registry_path', metavar='registry-path', help='path to the genome registry')
    init_parser.add_argument('--group-name',
        help='optional permission group name for write access to non-restricted subdirectories (e.g., user_defined_genes)')

    update_gene_parser = sp.add_parser('update-gene', help='update an existing user-defined gene with a new YAML gene model')
    update_gene_parser.set_defaults(func=update_user_defined_gene)
    update_gene_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    update_gene_parser.add_argument('--yaml-file', required=True,
        help = 'path to YAML file containing gene model for custom fasta')
    update_gene_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., hpc-durham, AWS, etc.)')

    get_genes_parser = sp.add_parser('get-genes',
        help='retrieve fasta and YAML genes models for selected user-defined genes')
    get_genes_parser.set_defaults(func=get_user_defined_genes)
    get_genes_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    get_genes_parser.add_argument('--gene-ids', required=True, nargs='+',
        help = 'comma-separated list of gene IDs to retrieve from registry')
    get_genes_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., hpc-durham, AWS, etc.)')
    get_genes_parser.add_argument('--outdir', required=False, default='./',
        help = 'target directory to write output to')
    get_genes_parser.add_argument('--version-delim', required=False, default='.',
        help = "delimiter when specifying a specific gene model version number after a gene-id (default = '.')")

    args = parser.parse_args()
    logger = start_logger(**vars(args))
    args.func(**vars(args))

    # TEST CODE
    # genome_list = []

    # base_genome = BaseGenome(
    #     id='mfas:106',
    #     version=106,
    #     species='mfas',
    #     assembly='Macaca_fascicularis_6.0',
    #     assembly_type='toplevel',
    #     sequence_type='dna',
    #     genome_fasta=GenomeFile(
    #         path={
    #             'hpc': '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/release-106/macaca_fascicularis/ensembl/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz'
    #         },
    #         type='fasta',
    #         source='genome'
    #     ),
    #     gtf=GenomeFile(
    #         path={
    #             'hpc': '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/release-106/macaca_fascicularis/ensembl/Macaca_fascicularis.Macaca_fascicularis_6.0.106.gtf.gz'
    #         },
    #         type='gtf'
    #     )
    # )
    # test_genome = Genome(
    #     id='mfas:106',
    #     base=base_genome,
    #     transcriptome_fasta=GenomeFile(
    #         path={
    #             'hpc': '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/release-106/macaca_fascicularis/derived/mfas.mfas60.ens106.transcriptome.toplevel.fa.gz'
    #         },
    #         type='fasta',
    #         source='transcriptome'
    #     ),
    #     star_index=GenomePath(
    #         path={
    #             'hpc': '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/release-106/macaca_fascicularis/star-index-271a'
    #         },
    #         type='star_index'
    #     ),
    #     refflat=GenomeFile(
    #         path={
    #             'hpc': '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/release-106/macaca_fascicularis/derived/mfas.mfas60.ens106.refflat'
    #         },
    #         type='refflat'
    #     ),
    #     rrna_interval_list=GenomeFile(
    #         path={
    #             'hpc': '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/release-106/macaca_fascicularis/derived/mfas.mfas60.ens106.rrna'
    #         },
    #         type='rrna_interval_list'
    #     )
    # )
    # genome_list.append(test_genome)

    # all_genomes = GenomeCollection(
    #     genomes={genome.id: genome for genome in genome_list}
    # )
 
    # conf_path = '/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/conf'
    # if not os.path.exists(conf_path):
    #     os.makedirs(conf_path)
    # with open(os.path.join(conf_path, '106.json'), 'w') as f:
    #     f.write(all_genomes.json())

    # with open(os.path.join(conf_path, '106.json'), 'r') as f:
    #     Genome.parse_file(f)
    #     model_json = json.load(f)
    # print(model_json)
    # model_json = Path(os.path.join(conf_path, '106.json'))
    # test_genome = Genome.parse_file(model_json)
    # print(test_genome.json())
    # new_genome = build_new_genome('/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/test.json')
    # print(new_genome.json())


    # register_genome(
    #     json_args='/lustre/workspace/home/moccir/rnaseq/nextflow/work/0b/f7a9e2e1bcd90eb6fcd619516628d7/test.json',
    #     conf_path='/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo')

    # register_genome(
    #     json_args='/lustre/workspace/home/moccir/rnaseq/nextflow/work/9d/c143f17df47aa5010e1fcdb52287ab/test.json',
    #     conf_path='/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo')

    # register_genome(
    #     json_args='/lustre/workspace/home/moccir/rnaseq/nextflow/work/4c/c3552f03ccc540a37fdf04ee17632e/test.json',
    #     conf_path='/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo')

    # register_genome(
    #     json_args='/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo/test.json',
    #     conf_path='/lustre/workspace/home/moccir/rnaseq/nextflow/test_repo')