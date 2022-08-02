import os
import pandas as pd

def get_selective_sweeps():
    with open ('highlight_0623.txt','r') as f1:
        list1=f1.readlines()
    with open('hs_3.gtf','r') as f2:
        list2=f2.readlines()
    with open('td_2.gtf','r') as f3:
        list3=f3.readlines()


    with open('selective_sweeps-td.txt','w') as f11:
        dict={}
        for i in list1:
            seq=[]
            i=i.strip()
            chr=i.split('\t')[0]
            start=i.split('\t')[1]
            end=i.split('\t')[2]

            for i4 in list3:
                i4=i4.strip()
                left=i4.split(' ')[2]
                right=i4.split(' ')[3]

                if chr in i4:
                    if int(left) > int(start) and int (right) < int(end):
                        id=i4.split(' ')[4]
                        seq.append(i)
                        if id not in dict:
                            dict[id] = seq
                        else:
                            dict[id] += seq

        for id,seq in dict.items():
            seq1 = [str(i) for i in seq]
            seq2 = list(set(seq1))
            seq2.sort(key=seq1.index)
            seq2 = ','.join(seq2)
            print(id + '\t'  + seq2, file=f11)

def merge():
     df1=pd.read_csv('hvat.txt',sep='\t')
     df2=pd.read_csv('selective_sweeps.txt',sep='\t')
     df3=pd.merge(df1,df2,how='inner',on='hvid')
     df3.to_csv('hvat_selective_sweeps.txt',sep='\t',index=None)

def main():
    os.chdir(']xpclr/')
    get_selective_sweeps()
    merge()
main()