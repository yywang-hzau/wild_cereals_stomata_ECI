import os
import pandas as pd

def get_big_block():
    with open('br.collinearity','r') as f1:
        with open('bigblock_bsbhs.txt','a+') as f2:
            with open('bigblock_bsbhd.txt','a+') as f3:
                list1=f1.readlines()
                list2=[]
                list3=[]
                for i in list1:
                    i=i.strip()
                    if '## Alignment' in i and  'BhD' in i and 'BhS' in i:
                        list2.append(i.split(' ')[2].split(':')[0]+'-')
                    if '## Alignment' in i and  'BhS' in i and 'BhD' in i:
                        list3.append(i.split(' ')[2].split(':')[0]+'-')
                for i1 in list1:
                    i1=i1.strip()
                    for i2 in list2:
                        if i2 in i1:
                            f2.write(i1+'\n')
                    for i3 in list3:
                        if i3 in i1:
                            f3.write(i1+'\n')

def rearrange_bh():
    f1=open ('br.gff','r')
    list1=f1.readlines()
    list2=[]
    list3=[]
    list4=[]
    list5=[]
    for i in list1:
        i=i.strip()
        if 'BhS' in i :
            list2.append(i.split('\t')[1])
        if 'BhD' in i:
            list3.append(i.split('\t')[1])
    f1.close()

    f2=open('bigblock_bsbhs.txt','r')
    for i2 in f2.readlines():
        list4.append(i2.strip().split('\t')[1])
        list4.append(i2.strip().split('\t')[2])
    list4 = set(list4)

    f3=open('bigblock_bsbhd.txt','r')
    for i3 in f3.readlines():
        list5.append(i3.strip().split('\t')[1])
        list5.append(i3.strip().split('\t')[2])
    list5=set(list5)

    f4=open('bhs_loss.txt','w')
    f5=open('bhd_loss.txt','w')

    for l2 in list2:
        if l2 not in list4:
            f4.write(l2+'\n')

    for l3 in list3:
        if l3 not in list5:
            f5.write(l3+'\n')

    f2.close()
    f3.close()
    f4.close()
    f5.close()

def additional_blast():
    f1=open('br.blast','r')
    list1=[]
    for i in f1.readlines():
        if 'Brahy.S' in i.split('\t')[0] and 'Brahy.D' in i.split('\t')[1]:
            list1.append(i.strip())
        if 'Brahy.S' in i.split('\t')[1] and 'Brahy.D' in i.split('\t')[0]:
            list1.append(i.strip())
    f2=open('bhs_loss.txt','r')
    f3=open('bhd_loss.txt','r')
    f4=open('bhs_add.txt','w')
    f5=open('bhd_add.txt','w')

    list2=f2.readlines()
    list3=f3.readlines()

    for i2 in list2:
        i2=i2.strip()
        for i1 in list1:
            if i2 in i1:
                f4.write(i1+'\n')
                continue
    for i3 in list3:
        i3=i3.strip()
        for i1 in list1:
            if i3 in i1:
                f5.write(i1+'\n')
                continue
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()

def cat():
    bhsbig=open('bigblock_bsbhs.txt','r')
    bhsadd=open('bhs_add.txt','r')
    bhs=open('bhd_blast.txt','w')
    bhs.write('block' + '\t' + 'bhdid' + '\t' + 'bhsid' + '\t' + 'evalue' + '\n')
    bhdbig=open('bigblock_bsbhd.txt','r')
    bhdadd=open('bhd_add.txt','r')
    bhd=open('bhs_blast.txt','w')
    bhd.write('block' + '\t' + 'bhsid' + '\t' + 'bhdid' + '\t' + 'evalue' + '\n')

    list1=bhsbig.readlines()
    list2=bhsadd.readlines()
    list3=bhdbig.readlines()
    list4=bhdadd.readlines()

    for i in list1:
        bhs.write(i)
    for i2 in list2:
        i2=i2.strip()
        n='nonblock'+'\t'+i2.split('\t')[0]+'\t'+i2.split('\t')[1]+'\t'+i2.split('\t')[10]
        bhs.write(n+'\n')
    for i in list3:
        bhs.write(i)
    for i2 in list4:
        i2=i2.strip()
        n='nonblock'+'\t'+i2.split('\t')[0]+'\t'+i2.split('\t')[1]+'\t'+i2.split('\t')[10]
        bhd.write(n+'\n')

    bhsbig.close()
    bhsadd.close()
    bhdbig.close()
    bhdadd.close()
    bhs.close()
    bhd.close()

def nohit():
    f1 = open('bhs_loss.txt', 'r')
    f2 = open('bhd_loss.txt', 'r')
    f3 = open('bhs_add.txt', 'r')
    f4 = open('bhd_add.txt', 'r')
    f5=open('bhs_nohit.txt','w')
    f6=open('bhd_nohit.txt','w')

    list1=f1.readlines()
    list2=f2.readlines()

    list3=[]
    for i3 in f3.readlines():
        list3.append(i3.strip().split('\t')[0])
    list3=set(list3)

    list4=[]
    for i4 in f4.readlines():
        list4.append(i4.strip().split('\t')[0])
    list4=set(list4)

    for i1 in list1:
        i1=i1.strip()
        if i1 not in list3:
            f5.write(i1 +'\n')

    for i2 in list2:
        i2=i2.strip()
        if i2 not in list4:
            f6.write(i2+'\n')

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()

def duplicate():
    df1= pd.read_csv('bhd_blast.txt', sep='\t')
    df3 = df1.sort_values(by='evalue', ascending=True)
    df5=df3.drop_duplicates(subset=['bhdid','bsid'],keep='first', inplace=False)
    df7 = df5.sort_values(by='block', ascending=True)
    bhdid = df7.drop_duplicates(['bhdid'], keep='first', inplace=False)

    df2=pd.read_csv('bigblock_bsbhs.txt',sep='\t')
    df4 = df2.sort_values(by='evalue', ascending=True)
    df6 = df4.drop_duplicates(subset=['bhsid', 'bhdid'], keep='first', inplace=False)
    df8 = df6.sort_values(by='block', ascending=True)
    bhdid=df8.drop_duplicates(['bhdid'], keep='first', inplace=False)
    bhsid = df8.drop_duplicates(['bhsid'], keep='first', inplace=False)

    df7.to_csv('bhd_drop.txt', sep='\t')
    df8.to_csv('bh_drop.txt', sep='\t')
    bhdid.to_csv('bhdunique.txt',sep='\t')
    bhsid.to_csv('bhsunique.txt', sep='\t')
    bhdid.to_csv('bhdunique.txt', sep='\t')

def value_count():
    f1=open('bh_drop.txt','r')
    f2=open('bhdunique.txt','r')
    f3=open('bhd_count.txt','a+')
    f3.write('bhdid'+'\t'+'block'+'\t'+'bhsid'+'\t'+'evalue'+'\n')
    list1=f1.readlines()
    list2=f2.readlines()
    list3=[]

    for i2 in list2[1:]:
        i2=i2.strip()
        list3.append(i2.split('\t')[2])

    for i3 in list3:
        count=0
        f3.write(i3+'\t')
        for i1 in list1:
            i1 = i1.strip()
            if i3 in i1:
                n=i1.split('\t')[1]+'\t'+i1.split('\t')[3]+'\t'+i1.split('\t')[4]
                f3.write(n+'\t'+'***'+'\t')
                count+=1
        f3.write(str(count)+'\n')


    f1.close()
    f2.close()
    f3.close()

    f4=open('bhd_drop.txt','r')
    f5=open('bhdunique.txt','r')
    f6=open('bhd_count.txt','a+')
    f6.write('bhid'+'\t'+'block'+'\t'+'bsid'+'\t'+'evalue'+'\n')
    list1=f4.readlines()
    list2=f5.readlines()
    list3=[]

    for i2 in list2[1:]:
        i2=i2.strip()
        list3.append(i2.split('\t')[2])

    for i3 in list3:
        count=0
        f6.write(i3+'\t')
        for i1 in list1:
            i1 = i1.strip()
            if i3 in i1:
                n=i1.split('\t')[1]+'\t'+i1.split('\t')[3]+'\t'+i1.split('\t')[4]
                f6.write(n+'\t'+'***'+'\t')
                count+=1
        f6.write(str(count)+'\n')

    f4.close()
    f5.close()
    f6.close()

def pick_block():
    f1=open('bhd_count.txt','r')
    f2=open('bhs_count.txt','r')
    list1=f1.readlines()
    list2=f2.readlines()
    f3=open('bhd_blockpick.txt','a+')
    f3.write('bhdid'+'\t'+'bsid'+'\t'+'add'+'\n')
    f4=open('bhs_blockpick.txt','a+')
    f4.write('bhsid' + '\t' + 'bsid' +'\t'+'add'+ '\n')

    for i1 in list1[1:]:
        i1=i1.strip()
        n=i1.split('\t')[-1]
        if 'nonblock' not in i1 and n =='1':
            i1new=i1.split('\t')[0]+'\t'+i1.split('\t')[2]+'\t'+i1.split('\t')[1]+i1.split('\t')[-1]+'\n'
            f3.write(i1new)

    for i2 in list2[1:]:
        i2=i2.strip()
        n = i2.split('\t')[-1]
        if 'nonblock' not in i2 and n == '1':
            i2new=i2.split('\t')[0]+'\t'+i2.split('\t')[2]+'\t'+i2.split('\t')[1]+i2.split('\t')[-1]+'\n'
            f4.write(i2new)

    f1.close()
    f2.close()
    f3.close()
    f4.close()

def merge_id():
    df1 = pd.read_csv('bhs_blockpick.txt', sep='\t')
    df2 = pd.read_csv('bhd_blockpick.txt', sep='\t')
    df3=pd.merge(df1,df2,on='bsid',how='inner')
    df3.to_csv('index.txt',sep='\t')

def get_fpkm_ct():
    df1 = pd.read_csv('bhaverage.txt', sep='\t')
    df1.insert(1, 'bhsid', df1['bhid'])
    df1.insert(2, 'bhdid', df1['bhid'])
    df2 = pd.read_csv('bsaverage.txt', sep='\t')
    column1 = df1.columns.values.tolist()
    column2 = df2.columns.values.tolist()
    bhc=[]
    bht=[]
    bsc=[]
    bst=[]
    for i1 in column1:
        if 'C' in i1 or 'id' in i1:
            bhc.append(i1)
        if 'T' in i1 or 'id' in i1:
            bht.append(i1)
    for i2 in column2:
        if 'C' in i2 or 'id' in i2:
            bsc.append(i2)
        if 'T' in i2 or 'id' in i2:
            bst.append(i2)
    df_bhc=df1[bhc]
    df_bht=df1[bht]
    df_bsc=df2[bsc]
    df_bst=df2[bst]
    df_bhc.to_csv('bhc.txt',sep='\t')
    df_bht.to_csv('bht.txt', sep='\t')
    df_bsc.to_csv('bsc.txt', sep='\t')
    df_bst.to_csv('bst.txt', sep='\t')

def merge_fpkm():
    df1 = pd.read_csv('bhc.txt', sep='\t')
    df2 = pd.read_csv('bht.txt', sep='\t')
    df3 = pd.read_csv('bsc.txt', sep='\t')
    df4 = pd.read_csv('bst.txt', sep='\t')
    df5 = pd.read_csv('index.txt', sep='\t')
    df6 = pd.merge(df1,df5,on='bhsid',how='inner')
    df7=pd.merge(df1,df5,on='bhdid',how='inner')
    df8=pd.merge(df3,df5,on='bsid',how='inner')
    df9 = pd.merge(df7, df8, on='bsid', how='inner')
    df10 = pd.merge(df6, df9, on='bsid', how='inner')
    df10.to_csv('cmerge.txt', sep='\t')

    df11 = pd.merge(df2, df5, on='bhsid', how='inner')
    df12= pd.merge(df2, df5, on='bhdid', how='inner')
    df13= pd.merge(df4, df5, on='bsid', how='inner')
    df14= pd.merge(df12, df13, on='bsid', how='inner')
    df15 = pd.merge(df11, df14, on='bsid', how='inner')
    df15.to_csv('tmerge.txt', sep='\t')

def trans_deg_metrix():
    df1=pd.read_csv('bsinner.csv', sep=',')
    df2=pd.read_csv('bhinner.csv', sep=',')
    df1.rename(columns={'xid': 'bsid'}, inplace=True)
    df2.rename(columns={'xid': 'bhid'}, inplace=True)
    df1=df1.fillna(0)
    df2=df2.fillna(0)
    column1 = df1.columns.values.tolist()
    column2 = df2.columns.values.tolist()

    for i1 in column1[1:]:
        df1[i1]=df1[i1].apply(lambda x: 0 if x==0 else 1)
    for i2 in column2[1:]:
        df2[i2]=df2[i2].apply(lambda x: 0 if x==0 else 1)
    df2.insert(1, 'bhsid', df2['bhid'])
    df2.insert(2, 'bhdid', df2['bhid'])
    df1.to_csv('bsdeg.txt', sep='\t')
    df2.to_csv('bhdeg.txt', sep='\t')



def main():
    os.chdir('/')
    #pick out big blocks from McScanX result
    get_big_block()
    #pick out genes for one of species that did not existed in big_block files
    rearrange_bh()
    #pick out additional blast hits for genes out of big blocks
    additional_blast()
    #combine blast result for each species
    cat()
    #pick out genes without optimum blast hits
    nohit()
    #drop duplicate in same block sorting by e-value
    duplicate()
    value_count()
    pick_block()
    merge_id()
    get_fpkm_ct()
    merge_fpkm()
    trans_deg_metrix()

main()
print('finish')