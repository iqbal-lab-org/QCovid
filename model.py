"""Declarative datamodel for covid assembly QC reporting pipeline
"""
from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()

class Dataset(Base):
    """Datasets map to ENA projects. Have many Runs.
    """
    __tablename__ = 'dataset'
    id = Column(Integer, primary_key=True)
    ena_id = Column(String(20))
    project_title = Column(String(200))
    runs = relationship("Run", back_populates="dataset")


class Run(Base):
    """A sequencing run for a single SAR-CoV-2 sample.
    """
    __tablename__ = 'run'
    id = Column(Integer, primary_key=True)
    ena_id = Column(String(20))
    dataset_id = Column(Integer, ForeignKey('dataset.id'))
    dataset = relationship("Dataset", uselist=False, back_populates="runs")

    assemblies = relationship("Assembly", back_populates="run")
    amplicon_qcs = relationship("AmpliconQC", back_populates="runs")

    dataset = relationship("Dataset", back_populates="runs")
    se_reads = relationship("SingleReads", uselist=False, back_populates="run")
    pe_reads = relationship("PairedReads", uselist=False, back_populates="run")


class SingleReads(Base):
    """Single end reads (i.e. Nanopore)
    """
    __tablename__ = 'se_reads'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    run = relationship("Run", back_populates="se_reads")
    read_length_mean = Column(Float)
    read_length_std = Column(Float)

    uri = Column(String(200))
    md5 = Column(String(32))


class PairedReads(Base):
    """Paired end reads (i.e. Illumina)
    """
    __tablename__ = 'pe_reads'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    run = relationship("Run", back_populates="pe_reads")

    read_count = Column(Integer)
    merged = Column(Integer)
    read_length_mean = Column(Float)
    read_length_std = Column(Float)

    r1_uri = Column(String(200))
    r1_md5 = Column(String(32))

    r2_uri = Column(String(200))
    r2_md5 = Column(String(32))


class Assembly(Base):
    """A complete 30k+ bp assembly derived from a sequencing run
    """
    __tablename__ = 'assembly'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    run = relationship("Run", back_populates="assemblies")

    name = Column(String(256))
    description = Column(String(500))
    sequence = Column(String(33000))

    mask = relationship("Mask", back_populates="assembly")


class SelfQC(Base):
    """Self consistency metrics for an assembly, given original reads
    """
    __tablename__ = 'self_qc'
    id = Column(Integer, primary_key=True)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))
    version = Column(String(16))
    f_95 = Column(Integer)
    f_90 = Column(Integer)
    f_85 = Column(Integer)
    f_80 = Column(Integer)
    f_75 = Column(Integer)
    f_00 = Column(Integer)
    ma_30 = Column(Integer)
    ns = Column(Integer)
    dashes = Column(Integer)
    unmapped = Column(Integer)
    bases = Column(Integer)
    template_size_mean = Column(Float)
    template_size_std = Column(Float)
    r1forward = Column(Integer)
    r2reverse = Column(Integer)
    secondary_alignments = Column(Integer)


class Mask(Base):
    """Actionable mask to be applied to poor regions of an assembly
    """
    __tablename__ = "mask"
    id = Column(Integer, primary_key=True)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))
    assembly = relationship("Assembly", back_populates="mask")
    criteria = Column(String(200))
    ref = Column(String(10))
    depth = Column(Integer)
    support = Column(Integer)
    position = Column(Integer)


class AmpliconQC(Base):
    """Performance of an individual amplicon in a sequencing run
    """
    __tablename__ = 'amplicon_qc'
    amplicon_id = Column(Integer, ForeignKey('amplicon.id'), primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'), primary_key=True)

    amplicons = relationship("Amplicon", back_populates="amplicon_qcs")
    runs = relationship("Run", back_populates="amplicon_qcs")

    reads = Column(Integer)
    forward_reads = Column(Integer)
    backward_reads = Column(Integer)
    unmapped_mates = Column(Integer)
    offtarget_mates = Column(Integer)
    reads_with_primer = Column(Integer)
    pairs_with_primers = Column(Integer)
    inner_fragments = Column(Integer)
    secondary_alignments = Column(Integer)
    read_crossing_interval = Column(Integer)
    fragment_crossing_interval = Column(Integer)
    version = Column(String(16))

class PrimerSet(Base):
    """Groups of amplicons
    """
    __tablename__ = 'primerset'
    id = Column(Integer, primary_key=True)
    name = Column(String(200))
    version = Column(String(20))
    reference = Column(String(100))

    amplicons = relationship("Amplicon", back_populates="primerset")


class Amplicon(Base):
    """Single amplicon from a primer set
    """
    __tablename__ = 'amplicon'
    id = Column(Integer, primary_key=True)
    primerset_id = Column(Integer, ForeignKey('primerset.id'))
    primerset = relationship("PrimerSet", back_populates="amplicons")
    amplicon_qcs = relationship("AmpliconQC", back_populates="amplicons")
    name = Column(String(50))
    fwd = Column(String(100))
    rev = Column(String(100))
    pool = Column(String(16))
    start = Column(Integer)
    end = Column(Integer)
