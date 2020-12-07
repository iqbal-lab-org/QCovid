import sys
from glob import glob
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from model import Dataset, Run, Assembly, PrimerSet, Amplicon, SelfQC, AmpliconQC, PairedReads, SingleReads
import model

VERSION = '0.0.1'

def illumina_metadata_batch1(session, p='../nov3/ILLUMINA_Metadata_Batch1.tsv'):
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

def nanopore_metadata(session, p='../oct5/OXFORD_NANOPORE_Metadata.tsv'):
    count = 0
    for l in open(p):
        project, sample, _, _, run_id, project_title, *rest = l.strip().split('\t')
        if project == 'project_id':
            continue

        prjs = list(session.query(Dataset).filter(Dataset.ena_id==project))
        if len(prjs) == 0:
            dataset = Dataset(ena_id=project, project_title=project_title)
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

    print(f"adding {count} nanopore runs")
    session.commit()


def nov3_assemblies(session, p='../nov3/*.fa'):
    count = 0
    for fasta in glob(p):
        description = None
        seq = None
        name = fasta.split('/')[-1].split('.')[0]
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
        run = session.query(Run).filter(Run.ena_id==sample).first()
        if run:
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
        name = fasta.split('/')[-1].split('.')[0]
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
        run = session.query(Run).filter(Run.ena_id==sample).first()
        if run:
            assembly = Assembly(description=description, name=name, sequence=seq, run_id=run.id)
            session.add(assembly)

            count += 1
    print(f"loaded {count} assemblies")
    session.commit()

def add_amplicons(session, bed, name=None):
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

def load_illumina_selfqc_batch1(session, p=''):
    for fn in glob.glob(p):
        pass

if __name__=="__main__":
    Session = sessionmaker()
    if len(sys.argv) == 3:
        if sys.argv[1] == 'init':
            engine = create_engine(sys.argv[2])

            Session.configure(bind=engine)
            model.Base.metadata.create_all(engine)
            session = Session()
            session.close()
            exit()
        else:
            print("unknown option")
            exit(1)
     
    elif len(sys.argv) < 2:
        engine = create_engine('sqlite://')

        model.Base.metadata.create_all(engine)
        Session.configure(bind=engine)

    else:
        engine = create_engine(sys.argv[1])
        Session.configure(bind=engine)
     
    session = Session()
    
#    nanopore_metadata(session)#, p='../oct5/smol_md.tsv')
#    illumina_metadata_batch1(session)#, p='../nov3/smol_md.tsv') 
   
#    add_amplicons(session, 'primers/nCoV-nl-primal500-75.bed')
#    add_amplicons(session, 'primers/nCoV-artic-v3.bed')

#    add_amplicon_qc(session, '../se_nl-primal.csv')
#    add_amplicon_qc(session, '../se_artic-v3.csv')
 #   oct5_assemblies(session)
 #   nov3_assemblies(session)

 #   add_self_qc(session, '../nl_cohort_qc.csv')

    add_pe(session)
    add_se(session)
    session.close()
