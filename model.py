from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Dataset(Base):
    __tablename__ = 'dataset'
    id = Column(Integer, primary_key=True)
    ena_id = Column(String(20))
    project_title = Column(String(200))
    runs = relationship("Run", back_populates="dataset")


class Run(Base):
    __tablename__ = 'run'
    id = Column(Integer, primary_key=True)
    ena_id = Column(String(20))
    dataset_id = Column(Integer, ForeignKey('dataset.id'))
    dataset = relationship("Dataset", uselist=False, back_populates="runs")

    assemblies = relationship("Assembly", back_populates="run")

    dataset = relationship("Dataset", back_populates="runs")
    se_reads = relationship("SingleReads", uselist=False, back_populates="run")
    pe_reads = relationship("PairedReads", uselist=False, back_populates="run")


class SingleReads(Base):
    __tablename__ = 'se_reads'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    run = relationship("Run", back_populates="se_reads")

    uri = Column(String(200))
    md5 = Column(String(32))


class PairedReads(Base):
    __tablename__ = 'pe_reads'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    run = relationship("Run", back_populates="pe_reads")

    r1_uri = Column(String(200))
    r1_md5 = Column(String(32))

    r2_uri = Column(String(200))
    r2_md5 = Column(String(32))


class Assembly(Base):
    __tablename__ = 'assembly'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    run = relationship("Run", back_populates="assemblies")

    name = Column(String(256))
    description = Column(String(500))
    sequence = Column(String(33000))

    mask = relationship("Mask", back_populates="assembly")


class SelfQC(Base):
    __tablename__ = 'self_qc'
    id = Column(Integer, primary_key=True)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))
    version = Column(String(16))
    f_95 = Column(Integer)
    f_90 = Column(Integer)
    #f_85 = Column(Integer)
    f_80 = Column(Integer)
    f_75 = Column(Integer)
    f_00 = Column(Integer)
    ns = Column(Integer)
    dashes = Column(Integer)
    bases = Column(Integer)


class Mask(Base):
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
    __tablename__ = 'amplicon_qc'
    amplicon_id = Column(Integer, ForeignKey('amplicon.id'), primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'), primary_key=True)
  
    amplicons = relationship("Amplicon", back_populates="amplicon_qcs")
    runs = relationship("Run", back_populates="amplicon_qcs")

    reads = Column(Integer)
    version = Column(String(16))


class PrimerSet(Base):
    __tablename__ = 'primerset'
    id = Column(Integer, primary_key=True)
    name = Column(String(200))
    version = Column(String(20))

    amplicons = relationship("Amplicon", back_populates="primerset")


class Amplicon(Base):
    __tablename__ = 'amplicon'
    id = Column(Integer, primary_key=True)
    primerset_id = Column(Integer, ForeignKey('primerset.id'))
    primerset = relationship("PrimerSet", back_populates="amplicons")
    name = Column(String(50))
    fwd = Column(String(100))
    rev = Column(String(100))
    pool = Column(String(16))
    start = Column(Integer)
    end = Column(Integer)
