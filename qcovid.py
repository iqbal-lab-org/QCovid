"""Query or populate QCovid database
"""
import sys
from os import path
from glob import glob

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from model import Dataset, Run, Assembly, PrimerSet, Amplicon, SelfQC, AmpliconQC, PairedReads, SingleReads
import model

VERSION = '0.0.1'

def crawl_ena_download(session, root, project):
    prj = session.query(Dataset).filter(Dataset.ena_id==project)
    prjroot = path.join(root, project)
    if not prj:
        print(f"Dataset {project} does not exist in database", file=sys.stderr)
        exit(1)
    if not path.isdir(prjroot):
        print(f"Dataset download directory does not exist: {prjroot}", file=sys.stderr)
        exit(1)

    # crawl single-read experiments downloaded into PRJ
    for sample in os.listdir(prjroot):
        # test that sample has not already been accessioned
        if '.txt' in exp:
            fq = os.listdir(path.join(prjroot, sample))
            if len(fq) == 1:
                # ensure _1.fq
                print(f"single,{project},{fq}")
                se += 1

            elif len(fq) == 2 or len(fq) == 3:
                fq1, fq2 = list(sorted(fq))[-2:]
                if '_1' not in fq1 or '_2' not in fq2:
                    # unexpected naming convention
                    continue
                pe += 1
                print(f"paired,{project},{fq1},{fq2}")

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

def do_se(session, out):
    no_reads = 0
    has_qc = 0
    score = 0
    for run in session.query(Run):
        se = session.query(SingleReads).filter(SingleReads.run_id==run.id).first()
        if not se:
            no_reads += 1
            continue

        for assembly in run.assemblies:
        #for assembly in session.query(Assembly).filter(Assembly.run_id==run.id):
            qc = session.query(SelfQC).filter(SelfQC.assembly_id==run.id).first()
            if qc:
                has_qc += 1
                continue
            score += 1
            dataset = session.query(Dataset).filter(Dataset.id==run.dataset_id).first()

            print(f"bsub.py 8 qcovid-self-{assembly.name} sh self_coverage_SE.sh {assembly.name} ../projects/{dataset.ena_id}/{run.ena_id}/{se.uri}", file=out)
        if (no_reads + has_qc + score) % 100 == 0:
            print(f"missing paired reads for {no_reads} runs. already have QC for {has_qc}. {score} PE reads with no QC")

parser = argparse.ArgumentParser(description='')
parser.add_argument('db')

subargs = parser.add_subparsers(dest='command')

init_args = subargs.add_parser('init')

prj_args = subargs.add_parser('project')

load_args = subargs.add_parser('load')
load_args.add_argument('dataset')

run_args = subargs.add_parser('run')

report_args = subargs.add_parser('report')

dump_args = subargs.add_parser('dump')

if __name__=="__main__":

    args = parser.parse_args()

    Session = sessionmaker()
    engine = create_engine(args.db)
    Session.configure(bind=engine)
     
    session = Session()

    if args.command == 'load':
        crawl_ena_download(session, "/tmp", sys.argv[2])
    elif parser.command == 'init':
        exit(1)
#    out = open('se_jobs.sh', 'w')
#    do_se(session, out) 
#    out.close() 
    session.close()
