from glob import glob
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from model import Dataset, Run, Assembly
import model

engine = create_engine('sqlite:///test.sqlite')

Session = sessionmaker()
Session.configure(bind=engine)
model.Base.metadata.create_all(engine)

session = Session()
md = {}
projects = set()
#for l in open('../nov3/ILLUMINA_Metadata_Batch1.tsv'):
for l in open('../nov3/smol_md.tsv'):
    _, project, sample, _, _, run, project_title, *rest = l.strip().split('\t')
    if project == 'sample_id':
        continue

    if project not in projects:
        dataset = Dataset(ena_id=project, project_title=project_title)
        session.add(dataset)
        projects.add(project)
    
    prj = list(session.query(Dataset).filter(Dataset.ena_id==project))
    assert len(prj_ids) == 1
    print(prj)
    run = Run(run_id=run, dataset_id=prj)
    prj.runs.add(run)
    session.add(run)

session.commit()

count = 0
for fasta in glob('../nov3/*.fa'):
    description = None
    seq = None
    name = fasta.split('/')[-1]
    sample = name.split('_')[0]
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

    rs = list(session.query(Run.id).filter(Run.run_id==name))
    if len(rs) == 0:
        continue
    assert len(rs) == 1


    sample = Run(ena_id=name)
    assembly = Assembly(description=description, seq=seq)
    session.add(assembly)

    count += 1
    if count > 30:
        break
session.commit()
session.close()
