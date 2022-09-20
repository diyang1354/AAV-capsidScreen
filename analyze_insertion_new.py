import subprocess
import pysam
import os
import gzip
import csv
import glob
import getopt
import sys
import re
import pandas as pd
from Bio import Seq
from Bio import SeqIO
def BBDukRemove(ref_fname, in_fname, out_fname, k,mkh):   
    #######filter out contaminants
    command = ['bbduk.sh',
               'ref=%s' % ref_fname,
               'in=%s' % in_fname,
               'out=%s' % out_fname,
               'k=%d' % k,
               'hdist=0',
               'mm=t', 
               'rcomp=t',
               'mkh=%d' % mkh
               ]
    ret = subprocess.call(command)
    assert ret == 0, 'bbduk.sh failed for %s!' % in_fname

def BBDukRemoveMulti(ref_fname, in_fnames, out_fnames, k,mkh):
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        BBDukRemove(ref_fname, in_fname, out_fname, k,mkh)

def BBDukRetain(ref_fname, in_fname, out_fname, k):   
    command = ['bbduk.sh',
               'ref=%s' % ref_fname,
               'in=%s' % in_fname,
               'outm=%s' % out_fname,
               'k=%d' % k,
               'hdist=0',
               'mm=t',
               'rcomp=f',
               'rename=t' ]
    ret = subprocess.call(command)
    assert ret == 0, 'bbduk.sh failed for %s!' % in_fname

def BBDukRetainMulti(ref_fname, in_fnames, out_fnames, k):
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        BBDukRetain(ref_fname, in_fname, out_fname, k)

def BBDukMask(ref_fname, in_fname, out_fname, k, kmask):
    command = ['bbduk.sh',
               'ref=%s' % ref_fname,
               'in=%s' % in_fname,
               'out=%s' % out_fname,
               'k=%d' % k,
               'kmask=%s' % kmask,
               'hdist=2',
               'mm=t',
               'rcomp=f']
    ret = subprocess.call(command)
    assert ret == 0, 'bbduk.sh failed for %s!' % in_fname

def BBDukMaskMulti(ref_fname, in_fnames, out_fnames, k, kmask):
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        BBDukMask(ref_fname, in_fname, out_fname, k, kmask)

def BBMapAlign(index_name, in_fname, out_fname):
    
    command = ['bbmap.sh',
               'ref=%s' % index_name,
               'in=%s' % in_fname,
               'out=%s' % out_fname,
               'touppercase=t',
               'local=f',
               'ambiguous=toss',
               'perfectmode=f',
               'nodisk']
    ret = subprocess.call(command)
    assert ret == 0, 'bbmap.sh failed for %s!' % in_fname

def SamtoolsIndex(in_fname):
    bam_fname = in_fname+'.bam'
    # Make a sorted BAM file.
    command = ['samtools view -bS %s' % in_fname,
               '|'
               'samtools sort -o %s' % bam_fname]
    ret = os.system(' '.join(command)) 
    assert ret == 0, 'samtools failed for %s!' % in_fname

    
    command = ['samtools index %s' % bam_fname]
    ret = os.system(' '.join(command))
    assert ret == 0, 'samtools failed to index %s!' % bam_fname

def calculate_insertion(aligned_reads_fnames):
    """ aligned_reads_fname : sorted bam file,
    just calculate foward insertion"""
    insertion_dic = {}
    #inspection_dic_left = {}
    #inspection_dic_right = {}
    for aligned_reads_fname in aligned_reads_fnames:
        aa_orient = re.search(r'.*\/.*_\d.(\S\S\S).\d.(.*).sam.bam',aligned_reads_fname,re.M|re.I).group(1)
        trim_orient = re.search(r'.*\/.*_\d.(\S\S\S).\d.(.*).sam.bam',aligned_reads_fname,re.M|re.I).group(2)
        alignedf = pysam.AlignmentFile(aligned_reads_fname, 'rb')
        for bbone_al in alignedf.fetch():
            backbone_match_strand = -1 if bbone_al.is_reverse else 1
            if(aa_orient == 'fwd'  and backbone_match_strand >0):
                if(trim_orient == 'left_trimmed' and bbone_al.reference_start != 0):  
                    insertion_index = bbone_al.reference_start+5
                    if insertion_dic.__contains__(str(insertion_index)):
                        insertion_dic[str(insertion_index)] += 1
                    else:
                        insertion_dic[str(insertion_index)] = 1
                elif(trim_orient == 'right_trimmed' and bbone_al.reference_end != 1604):
                    insertion_index = bbone_al.reference_end
                    if insertion_dic.__contains__(str(insertion_index)):
                        insertion_dic[str(insertion_index)] += 1
                    else:
                        insertion_dic[str(insertion_index)] = 1   
                else:
                    continue            
            elif(aa_orient == 'rev'  and backbone_match_strand <0):
                if(trim_orient == 'left_trimmed' and bbone_al.reference_end != 1604):  
                    insertion_index = bbone_al.reference_end
                    if insertion_dic.__contains__(str(insertion_index)):
                        insertion_dic[str(insertion_index)] += 1
                    else:
                        insertion_dic[str(insertion_index)] = 1
                elif(trim_orient == 'right_trimmed' and bbone_al.reference_start != 0 ):
                    insertion_index = bbone_al.reference_start+5
                    #inspection_dic_right[bbone_al.qname] = bbone_al                    
                    if insertion_dic.__contains__(str(insertion_index)):
                        insertion_dic[str(insertion_index)] += 1
                    else:
                        insertion_dic[str(insertion_index)] = 1
                else:
                    continue
            else:
                continue
    filter_dic={}
    for key in insertion_dic.keys():
        if int(key)%3 == 0:
            filter_dic[str(int(int(key)/3+202))] = insertion_dic[key]
    return filter_dic
def main(argv):
    ''' read in arguments '''
    try:
        opts, args = getopt.gnu_getopt(argv,"",['help','read_1=','read_2=','insert_fname=','backbone_fname=','tmp_file=','out=','insert_len=','kmer_remove=','kmer_retain='])
    except getopt.GetoptError:
        print("option not found! Exit now")
        sys.exit(2)
        
    for opt_name,opt_value in opts:
        if opt_name == "--help":
            print("Analyze insertion positions and their counts\nParameters:\n--read_1:\tpath to read_1\n--read_2:\tpath to read_2\n--backbone_fname:\tpath to backbone reference\n--insert_fname:\tpath to insert reference\n--tmp_file:\tpath to save temporary files\n--out:\tpath to output csv file containning positions and counts,only write path without file itself\n--insert_len:\tinsertion length\n--kmer_remove:\tkmers to remove multi-insert\n--kmer_ratain:\tkmers to retain insertion")
            exit()
        elif opt_name == '--read_1':
            fastq_R1 = opt_value
            
        elif opt_name == '--read_2':
            fastq_R2 = opt_value
            
        elif opt_name == '--insert_fname':
            ref_fname = opt_value
            
        elif opt_name == '--tmp_file':
            tmp_path = opt_value
        
        elif opt_name == '--backbone_fname':
            bbone_ref = opt_value

        elif opt_name == '--insert_len':
            insert_ref_len = int(opt_value)
        
        elif opt_name == '--kmer_remove':
            kmer_remove = int(opt_value)
        
        elif opt_name == '--kmer_retain':
            kmer_retain = int(opt_value)

        elif opt_name == '--out':
            out_csv_path = opt_value
    base1=os.path.basename(fastq_R1).split('.')[0]
    base2=os.path.basename(fastq_R2).split('.')[0]
    """ filter multiple-inserts reads"""
    #insert_ref_len = 37
    #kmer_retain = 20
    #kmer_remove = 12
    min_kmer_hit = insert_ref_len - kmer_remove + 2  #### used to remove multiple insert
    print('filtering  multiple inserts  \n')
    multi_insert_filtered_R1=os.path.join(tmp_path, base1+'.multifiltered.fq.gz')
    multi_insert_filtered_R2=os.path.join(tmp_path, base2+'.multifiltered.fq.gz')
    BBDukRemoveMulti(ref_fname=ref_fname, in_fnames=[fastq_R1,fastq_R2], out_fnames=[multi_insert_filtered_R1,multi_insert_filtered_R2],k=kmer_remove,mkh=min_kmer_hit) ### 

    print("done! \n Now filter reads contain insert")
    insert_filtered_R1 = os.path.join(tmp_path, base1+'.aaFiltered.fq.gz')
    insert_filtered_R2 = os.path.join(tmp_path, base2+'.aaFiltered.fq.gz')
    
    BBDukRetainMulti(ref_fname=ref_fname,in_fnames=[multi_insert_filtered_R1,multi_insert_filtered_R2],out_fnames=[insert_filtered_R1,insert_filtered_R2],k=kmer_retain)
    print("done! \n")

    print("Mask 5' end of insert in the reads \n")
    masked_R1 = os.path.join(tmp_path, base1+'.masked.fq.gz')
    masked_R2 = os.path.join(tmp_path, base2+'.masked.fq.gz')
    BBDukMaskMulti(ref_fname=ref_fname,in_fnames=[insert_filtered_R1,insert_filtered_R2],out_fnames=[masked_R1,masked_R2],k=kmer_retain,kmask='Z')
    print("done! \n Start to trim file")

    for file_name in [masked_R1,masked_R2]:
        base = os.path.basename(file_name).split('.')[0]
        os.system(' '.join(["zcat",file_name ,"|", "grep",'rev', "|", "sed", 's/@//g' ,">", os.path.join(tmp_path,base+".rev.txt")]))
        os.system(' '.join(["zcat",file_name ,"|", "grep",'fwd', "|", "sed", 's/@//g' ,">", os.path.join(tmp_path,base+".fwd.txt")]))
        masked_fwd = os.path.join(tmp_path,base+'.fwd.fq.gz')
        masked_rev = os.path.join(tmp_path,base+'.rev.fq.gz')
        os.system(' '.join(["filterbyname.sh","in=%s" %file_name, "out=%s" % masked_fwd, "names=%s" % os.path.join(tmp_path,base+".fwd.txt"), "include=t","ignorejunk=t"]))
        os.system(' '.join(["filterbyname.sh","in=%s" %file_name, "out=%s" % masked_rev, "names=%s" % os.path.join(tmp_path,base+".rev.txt"), "include=t","ignorejunk=t"]))
        for masked_file in [masked_fwd,masked_rev]:
            prefix = re.search(r'.*\/(.*).fq.gz',masked_file).group(1)
            print('Start to trim file: %s' %prefix)
            left_trimmed = []
            right_trimmed = [] 
            with gzip.open(masked_file,"rt") as handle:
                for record in SeqIO.parse(handle,'fastq'):
                    z_pos_first = str(record.seq).find('Z')
                    if z_pos_first >= 20:
                        record= record[0:z_pos_first]
                        right_trimmed.append(record)
                    else:
                        z_pos_last = str(record.seq).rfind('Z')
                        record = record[z_pos_last+1:]
                        left_trimmed.append(record)
            SeqIO.write(left_trimmed,os.path.join(tmp_path,prefix+'.left_trimmed.fq'),"fastq")
            SeqIO.write(right_trimmed,os.path.join(tmp_path,prefix+'.right_trimmed.fq'),"fastq")
            print('Finished left and right trimming: %s' %prefix)

    for file_name in glob.glob(os.path.join(tmp_path,'*trimmed.fq')):
        base_name = os.path.basename(file_name).split('.')[0]
        aa_orient = re.search(r'.*\/.*_(\d).(\S\S\S).(.*).fq',file_name,re.M|re.I).group(2)
        read_num = re.search(r'.*\/.*_(\d).(\S\S\S).(.*).fq',file_name,re.M|re.I).group(1)
        trim_orient = re.search(r'.*\/.*_(\d).(\S\S\S).(.*).fq',file_name,re.M|re.I).group(3)
        alignment = os.path.join(tmp_path,base_name+'.'+aa_orient+'.'+read_num+'.'+trim_orient+'.sam')
        BBMapAlign(index_name=bbone_ref,in_fname= file_name, out_fname=alignment)
        SamtoolsIndex(alignment)
    R1_files= glob.glob(os.path.join(tmp_path,'*.*.*.*.bam'))
    #R2_files= glob.glob(os.path.join(tmp_path,'*.*.2.*.bam'))
    diction_R1= calculate_insertion(R1_files)
    #diction_R2= calculate_insertion(R2_files)
    file_base = base_name.split('_')[0]
    df = pd.DataFrame(data=diction_R1, index=[file_base])
    df = (df.T)
    df.to_excel(os.path.join(out_csv_path,file_base+'_results.xlsx'))
    os.system(' '.join(['rm',tmp_path+'/'+'*']))

if __name__ == "__main__":
    main(sys.argv[1:])
   
