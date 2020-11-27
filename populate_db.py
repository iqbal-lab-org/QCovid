import sys
from glob import glob
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from model import Dataset, Run, Assembly
import model


def illumina_metadata_batch1(session, p='../nov3/ILLUMINA_Metadata_Batch1.tsv'):
    projects = set()
    count = 0
    for l in open(p):
        _, project, sample, _, _, run, project_title, *rest = l.strip().split('\t')
        if project == 'sample_id':
            continue

        if project not in projects:
            dataset = Dataset(ena_id=project, project_title=project_title)
            session.add(dataset)
            projects.add(project)
            session.commit()
        
        prjs = list(session.query(Dataset).filter(Dataset.ena_id==project))
        assert len(prjs) == 1
        prj = prjs[0]
        run = Run(run_id=run, dataset_id=prj)
        prj.runs.append(run)
        count += 1
        session.add(run)
    print(f"adding {count} illumina runs")
    session.commit()

def nanopore_metadata(session, p='../oct5/OXFORD_NANOPORE_Metadata.tsv'):
    projects = set()
    count = 0
    for l in open(p):
        project, sample, _, _, run_id, project_title, *rest = l.strip().split('\t')
        if project == 'project_id':
            continue

        if project not in projects:
            dataset = Dataset(ena_id=project, project_title=project_title)
            session.add(dataset)
            projects.add(project)
            session.commit()
        
        prjs = list(session.query(Dataset).filter(Dataset.ena_id==project))
        if len(prjs) != 1:
            continue
        prj = prjs[0]
        run = Run(run_id=run_id, dataset_id=prj)
        prj.runs.append(run)
        count += 1
        session.add(run)
    print(f"adding {count} nanopore runs")
    session.commit()


def nov3_assemblies(session, p='../nov3/*.fa'):
    count = 0
    for fasta in glob(p):
        description = None
        seq = None
        name = fasta.split('/')[-1]
        sample = name.split('.')[0].split('_')[0]
        with open(fasta) as fd:
            ls = 0
            for l in fd:
                ls += 1
                if l[0] == '>':
                    description = l.rstrip()
                else:
                    seq = l.rstrip()
            if ls > 2:
                print("error")
                continue
        rs = list(session.query(Run).filter(Run.run_id==sample))
        if len(rs) == 0:
            continue
        assert len(rs) == 1

        run = rs[0]
        #sample = Run(ena_id=name)
        assembly = Assembly(description=description, name=name, sequence=seq, run_id=run.id)
        session.add(assembly)

        count += 1
    print(f"loaded {count} assemblies")
    session.commit()

def oct5_assemblies(session, p='../oct5/*.fasta'):
    count = 0
    for fasta in glob(p):
        description = None
        seq = None
        name = fasta.split('/')[-1]
        sample = name.split('.')[0].split('_')[0]
        with open(fasta) as fd:
            ls = 0
            for l in fd:
                ls += 1
                if l[0] == '>':
                    description = l.rstrip()
                else:
                    seq = l.rstrip()
            if ls > 2:
                print("error")
                continue
        rs = list(session.query(Run).filter(Run.run_id==sample))
        if len(rs) == 0:
            continue
        assert len(rs) == 1

        run = rs[0]
        #sample = Run(ena_id=name)
        assembly = Assembly(description=description, name=name, sequence=seq, run_id=run.id)
        session.add(assembly)

        count += 1
    print(f"loaded {count} assemblies")
    session.commit()


def load_illumina_selfqc_batch1(session, p=''):
    for fn in glob.glob(p):
        pass

if __name__=="__main__":
    if len(sys.argv) < 2:
        db = 'sqlite://'
    else:
        db = sys.argv[1]
    engine = create_engine(db)

    Session = sessionmaker()
    Session.configure(bind=engine)
    model.Base.metadata.create_all(engine)
    
    session = Session()
    
#    illumina_metadata_batch1(session)#, p='../nov3/smol_md.tsv') 
#    nanopore_metadata(session)
    
#    oct5_assemblies(session)

#    nov3_assemblies(session)
    session.close()
