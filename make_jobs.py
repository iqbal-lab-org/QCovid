import sys
from glob import glob
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from model import Dataset, Run, Assembly, PrimerSet, Amplicon, SelfQC, AmpliconQC, PairedReads, SingleReads
import model

VERSION = '0.0.1'

def dump_assembly(session, aid):
    fa = open(f"{assembly.name}.fa", 'w')
    print(f">{assembly.description}\n{assembly.sequence}", file=fa)
    fa.close()
    return f"{assembly.name}.fa"

def do_pe(session, out):
    no_reads = 0
    has_qc = 0
    score = 0
    for run in session.query(Run):
        pe = session.query(PairedReads).filter(PairedReads.run_id==run.id).first()
        if not pe:
            no_reads += 1
            continue

        for assembly in session.query(Assembly).filter(Assembly.run_id==run.id):
            qc = session.query(SelfQC).filter(SelfQC.assembly_id==run.id).first()
            if qc:
                has_qc += 1
                continue
            score += 1
            dataset = session.query(Dataset).filter(Dataset.id==run.dataset_id).first()

            print(f"bsub.py 8 qcovid-self-{assembly.name} sh self_coverage_PE.sh {assembly.name} ../projects/{dataset.ena_id}/{run.ena_id}/{pe.r1_uri} ../projects/{dataset.ena_id}/{run.ena_id}/{pe.r2_uri}", file=out)
        if (no_reads + has_qc + score) % 100 == 0:
            print(f"missing paired reads for {no_reads} runs. already have QC for {has_qc}. {score} PE reads with no QC")

if __name__=="__main__":
    Session = sessionmaker()
    engine = create_engine(sys.argv[1])
    Session.configure(bind=engine)
     
    session = Session()
   
    out = open('pe_jobs.sh', 'w')
    do_pe(session, out) 
    out.close() 
    session.close()
