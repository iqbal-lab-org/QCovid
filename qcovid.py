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

VERSION = '0.0.2'

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

def dump_assembly(session):
    fa = open(f"{assembly.name}.fa", 'w')
    print(f">{assembly.description}\n{assembly.sequence}", file=fa)
    fa.close()
    return f"{assembly.name}.fa"

def run_pe(session, out):
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

def run_se(session, out):
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

def import_amplicons(session, bed, name=None):
    """Import amplicons and metadata for a primer set
    """
    amplicons = []
    if name is None:
        name = bed.split('.')[0].split('/')[-1]

    if session.query(PrimerSet).filter(PrimerSet.name==name).first():
        print(f"primerset {name} already loaded. skipping.")
        return

    pset = PrimerSet(reference='MN908947', name=name, version=VERSION)
    session.add(pset)
    session.commit()
    count = 0
    for line in open(bed):
        if line[0] == '#':
            continue

        n, start, end = line.strip().split('\t')
        start = int(start)
        end = int(end)
        amplicon = Amplicon(primerset_id=pset.id, name=n, start=start, end=end)
        session.add(amplicon)
        count += 1

    print(f"loaded primerset {name}, {count} amplicons")
    session.commit()


def import_metadata_batch(session, p):
    """Historical: load project metadata from a TSV file
    """
    count = 0
    for l in open(p):
        _, project, sample, _, _, run_id, project_title, *rest = l.strip().split('\t')
        if project == 'sample_id':
            continue
        
        prjs = list(session.query(Dataset).filter(Dataset.ena_id==project))
        if len(prjs) == 0:
            dataset = Dataset(ena_id=project, project_title=project_title)
            prj = dataset
            session.add(dataset)
            session.commit()
        prjs = list(session.query(Dataset).filter(Dataset.ena_id==project))

        assert len(prjs) == 1
        prj = prjs[0]
       
        runq = session.query(Run).filter_by(ena_id=run_id).first()
        if not runq:
            run = Run(ena_id=run_id, dataset_id=prj)
            prj.runs.append(run)
            count += 1
            session.add(run)

    print(f"adding {count} illumina runs")
    session.commit()


parser = argparse.ArgumentParser(description='QCovid database interface tools')
parser.add_argument('db')

subargs = parser.add_subparsers(dest='command')

init_args = subargs.add_parser('init')

prj_args = subargs.add_parser('project')

load_args = subargs.add_parser('load')
load_args.add_argument('dataset')
load_args.add_argument('-d', default='./')

run_args = subargs.add_parser('run')

report_args = subargs.add_parser('report')

dump_args = subargs.add_parser('dump')

info_args = subargs.add_parser('info')

if __name__=="__main__":

    args = parser.parse_args()

    Session = sessionmaker()
    engine = create_engine(args.db)
    Session.configure(bind=engine)

    if args.command == 'load':
        """Crawl ena_data_get directory and import fastq file pairs into database"""
        session = Session()
        crawl_ena_download(session, load_args.d, load_args.dataset)

    elif parser.command == 'init':
        """Initialise a QCovid database instance"""
        model.Base.metadata.create_all(engine)
        session = Session()

    elif parser.command == 'import':
        """historical: load data from raw tsv file"""
        import_metadata_batch(session, p=import_args.tsv)

    elif parser.command == 'run':
        """For each sample with fastq data but no QC data, run pipeline""" 
        pass

    elif parser.command == 'info':
        """Dump basic stats from database connection"""

    session.close()


def add_self_qc(session, fn):
    count = 0
    for line in open(fn):
        if line[0] == '#':
            continue
        sample, f_95, f_90, f_80, f_75, f_00, ma_30, ns, dashes = line.strip().split(',')

        run = session.query(Run).filter(Run.ena_id==sample).first()
        if not run:
            continue
        assembly = session.query(Assembly).filter(Assembly.run_id==run.id).first()
        if assembly:
            qcs = SelfQC(assembly_id=assembly.id, version=VERSION, f_95=int(f_95),\
                    f_90=int(f_90), f_80=int(f_80), f_75=int(f_75), f_00=int(f_00),\
                    ma_30 = int(ma_30), ns=int(ns), dashes=int(dashes))
            session.add(qcs)
            count += 1
    print(f"added {count} lines of self-QC")
    session.commit()

def add_pe(session, fn='../crawl_prj_PE.o', pes='../nov3/ILLUMINA_Metadata_Batch1.tsv'):
    runs = set()
    for l in open(pes):
        runs.add(l.strip().split('\t')[5])

    count = 0
    for line in open(fn):
        prj, run, fq1, fq2 = line.strip().split(',')
        if run not in runs:
            continue

        run = session.query(Run).filter(Run.ena_id==run).first()
        if not run:
            continue

        pe_ = session.query(PairedReads).filter(PairedReads.run_id==run.id).first()
        if pe_:
            continue

        pe = PairedReads(run_id=run.id, r1_uri=fq1, r2_uri=fq2)
        session.add(pe)
        count += 1
    print(f"added {count} paired end reads")
    session.commit()

def add_se(session, fn='../crawl_prj_SE.o', ses='../oct5/OXFORD_NANOPORE_Metadata.tsv'):
    runs = set()
    for l in open(ses):
        runs.add(l.strip().split('\t')[4])

    count = 0
    for line in open(fn):
        prj, run, fq = line.strip().split(',')
        fq = fq[2:-2]
        if run not in runs:
            continue

        run = session.query(Run).filter(Run.ena_id==run).first()
        if not run:
            continue

        se_ = session.query(SingleReads).filter(SingleReads.run_id==run.id).first()
        if se_:
            continue

        se = SingleReads(run_id=run.id, uri=fq)
        session.add(se)
        count += 1
    print(f"added {count} single end reads")
    session.commit()


def add_amplicon_qc(session, fn, primerset=None):
    count = 0
    skipped = 0
    header = None
    for line in open(fn):
        if line[0] == '#':
            _, *header = line.strip().split(',')
            continue

        sample, *depths = line.strip().split(',')

        run = session.query(Run).filter(Run.ena_id==sample).first()
        if not run:
            continue

        for amplicon_name, depth in zip(header, depths):
            depth = int(depth)
            amplicon = session.query(Amplicon).filter(Amplicon.name==amplicon_name).first()
            assert amplicon
            if session.query(AmpliconQC).filter(AmpliconQC.amplicon_id==amplicon.id, AmpliconQC.run_id==run.id, AmpliconQC.version==VERSION).first():
                skipped += 1
                continue
            aqc = AmpliconQC(amplicon_id=amplicon.id, run_id=run.id, reads=depth, version=VERSION)
            session.add(aqc)
        session.commit()
        count += 1
    print(f"added amplicon QC for {count} runs (skipped {skipped} QC rows)")
