from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Dataset(Base):
    __tablename__ = 'dataset'
    id = Column(Integer, primary_key=True)
    ena_id = Column(String)
    project_title = Column(String)
    runs = relationship("Run", back_populates="dataset")


class Platform(Base):
    __tablename__ = 'platform'
    id = Column(Integer, primary_key=True)


class Run(Base):
    __tablename__ = 'run'
    id = Column(Integer, primary_key=True)
    run_id = Column(String)
    dataset_id = Column(Integer, ForeignKey('dataset.id'))
    #assemblies = relationship("Assembly", back_populates="run")
    dataset = relationship("Dataset", back_populates="runs")

class SingleReads(Base):
    __tablename__ = 'se_reads'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    assembly_id = Column(Integer, ForeignKey('assembly.id'))

    uri = Column(String)
    md5 = Column(String)


class PairedReads(Base):
    __tablename__ = 'pe_reads'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    assembly_id = Column(Integer, ForeignKey('assembly.id'))

    r1_uri = Column(String)
    r1_md5 = Column(String)

    r2_uri = Column(String)
    r2_md5 = Column(String)


class Assembly(Base):
    __tablename__ = 'assembly'
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    name = Column(String)
    description = Column(String)
    sequence = Column(String)

class SelfQC(Base):
    __tablename__ = 'self_qc'
    id = Column(Integer, primary_key=True)
    assembly_id = Column(Integer, ForeignKey('assembly.id'))


class AmpliconQC(Base):
    __tablename__ = 'amplicon_qc'
    id = Column(Integer, primary_key=True)
    

class PrimerSet(Base):
    __tablename__ = 'primers'
    id = Column(Integer, primary_key=True)
    name = Column(String)
