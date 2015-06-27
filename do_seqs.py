# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:14:44 2015

@author: Tiago
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
Entrez.email = "tiago_alves26@hotmail.com"


class do_seqs:

    def __init__(self,term=0,mutations={}):
        self.term=term #usado para criar o ficheiro quando o input e uma proteina
        self.mutations=mutations #posicoes e alteracao que ocorre aquando de uma mutacao        
        self.seq="" #sequencia normal
        self.mseq="" #sequencia mutada
    
       
###############################################################################
       #Pesquisa por gene e posterior conversao em proteina
###############################################################################
       
       
    def get_gene_file(self,gene=0,idn=0):
        '''
        Funcao que dado um GI ou um nome de um gene, 
        extrai a sequencia da base dados no NCBI e guarda-a num ficheiro FASTA, retornando o seu nome.
        '''
        return_filename=""        
        if idn==0:#damos o nome do gene a funcao 
            hand=Entrez.esearch(db='gene',term=gene+"[sym]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            
            handl=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=idnum)
            record=Entrez.read(handl)
            if 'Entrezgene_comments' in record[0].keys():
                reg=record[0]['Entrezgene_comments']
                for pos in range(len(reg)):
                    if 'Gene-commentary_comment' in record[0]['Entrezgene_comments'][pos].keys():
                        reg2=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0].keys()
                        if 'Gene-commentary_comment' in reg2:
                            reg3=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]
                            for i in range(len(reg3)):
                                try:
                                    to=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para
                                    desde=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de
                                    identif=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                                except:
                                    pass
            
                            if 'Gene-commentary_comment' in reg3:
                                reg4=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment']
                                for i in range(len(reg4)):
                                    try:
                                        to=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para
                                        desde=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de
                                        identif=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                                    except:
                                        pass
                        else:
                            pass
                                
                    else:
                        pass
                    
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )#depois de termos o id, e os valores onde comeca e acaba a seq vamos busca-la
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(gene).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
          
        elif gene==0:#damos o gi do gene
            hand=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=str(idn))
            record=Entrez.read(hand)
            if 'Entrezgene_comments' in record[0].keys():
                reg=record[0]['Entrezgene_comments']
                for pos in range(len(reg)):
                    if 'Gene-commentary_comment' in record[0]['Entrezgene_comments'][pos].keys():
                        reg2=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0].keys()
                        if 'Gene-commentary_comment' in reg2:
                            reg3=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]
                            for i in range(len(reg3)):
                                try:
                                    to=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para
                                    desde=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de
                                    identif=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                                except:
                                    pass
            
                            if 'Gene-commentary_comment' in reg3:
                                reg4=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment']
                                for i in range(len(reg4)):
                                    try:
                                        to=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para
                                        desde=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de
                                        identif=record[0]['Entrezgene_comments'][pos]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_comment'][i]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                                    except:
                                        pass
                        else:
                            pass
                                
                    else:
                        pass
                    
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )#depois de termos o id, e os valores onde comeca e acaba a seq vamos busca-la
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(idn).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
        
        return return_filename


                        
    def get_seq_from_gene(self,gene=0,idn=0):
        '''
        Funcao que dado um GI ou um nome de um gene, 
        vai retornar a sequencia desse dado input, em formato de texto, criando um ficheiro chamando a funcao get_gene_file
        '''
        if idn==0:
            gf=self.get_gene_file(gene)
            self.term=gene #atribuir o termo ao pesquisado para quando converter em proteina saber qual o ficheiro a abrir
        elif gene==0:
            gf=self.get_gene_file(0,idn)
            self.term=idn
        f=open(str(gf), 'r')
        sequence = SeqIO.read(f,"fasta",generic_dna)
        self.seq=str(sequence.seq) #atualizar a classe para esta sequencia
        print ("A sequencia do gene pesquisado é:\n"+str(sequence.id)) #para informar qual sera a sequencia usada
        return str(sequence.seq)
    
    
    
    def convert_into_protein(self,seq):#so se chama depois das anteriores para ter o termo associado a pesquisa efetuada
        '''
        Funcao que dada uma sequencia de DNA, a converte numa proteína e cria um ficheiro com a mesma em formato FASTA.
        '''        
        s=Seq(seq,generic_dna)
        trans=s.translate()
        self.seq=trans#sequencia normal
        fil=open('protein_'+str(self.term)+'.fasta','w+')#grava novo ficheiro
        fil.writelines(trans)
        fil.close()
        return trans
        
    
    
###############################################################################
       #Pesquisa nas bases de dados sabendo qual a proteina a pesquisar
###############################################################################
    
    def search_seq_ids(self,term,db="protein"):
        '''
        Funcao procura GIs dado um termo a pesquisar na base de dados 'Protein'
        '''
        handle=Entrez.esearch(db,term)
        record=Entrez.read(handle)
        ids=record["IdList"]
        return ids
        
        
    def get_seq(self,id_number=0,term=0):
        '''
        Funcao que recorre ao NCBI para retirar a sequencia de uma proteina, 
        dado o seu GI ou o termo associado a ela
        '''
        if id_number==0:
            handle=Entrez.esearch('protein',term)
            record=Entrez.read(handle)
            idnum=record["IdList"][0]#primeiro elemento da lista de resultados
            handl=Entrez.efetch(db='protein',rettype="gb",retmode="xml",id=idnum)
            results=Entrez.read(handl)
            sequence=results[0]["GBSeq_sequence"]#sequencia
            self.seq=sequence.upper()
            definition=results[0]['GBSeq_definition']
            name=term
            print ("Searched protein sequence belongs to:\n"+str(definition)+'\nGI:'+idnum)#para sabermos de que proteina se trata        
        elif term==0:
            handle = Entrez.efetch(db="protein", id=id_number, rettype="gb", retmode="xml")
            results=Entrez.read(handle)
            sequence=results[0]["GBSeq_sequence"]
            self.seq=sequence.upper()        
            definition=results[0]['GBSeq_definition']
            name=id_number
            print ("Searched protein sequence belongs to:\n"+str(definition)+'\nGI:'+id_number)
       
        return self.seq,name
     
     
     
##corrida normal, dada a proteina.     
    def create_file_with_seq(self,id_number=0,term=0,mutated=False):
        '''
        Funcao que cria ficheiros com a proteina (estado normal e estado mutado), 
        sendo que permite a pesquisa de proteinas no estado normal (pelo GI ou termo)
        
        '''
        if id_number==0:
            if not mutated:    
                seq,name=self.get_seq(0,term)
                fich=open('protein_'+str(name)+'.fasta','w+')
                fich.write(seq)
            else:
                seq=self.mseq
                fich=open('protein_'+str(term)+'mutated.fasta','w+')
                fich.write(seq)
        
        elif term==0:
            if not mutated:    
                seq,name=self.get_seq(id_number)
                fich=open('protein_'+str(name)+'.fasta','w+')
                fich.write(seq)
            else:
                seq=self.mseq
                fich=open('protein_'+str(id_number)+'mutated.fasta','w+')                
                fich.write(seq)
        
        fich.close()


            
#corrida mutada, apos , chamar a anterior para criar o ficheiro
            
    def induce_mutation(self):
        '''
        Funcao que insere as mutacoes na sequencia normal
        '''
        seq_m=""
        for pos in self.mutations.keys():
            if self.seq[int(pos)-1]==self.mutations[pos][0]:#primeiro elemento do tuplo (am. inicial, am. final)
                seq_m+=self.seq[:int(pos)-1]
                seq_m+=str(self.mutations[pos][1])#aminoacido de substituicao
                seq_m+=self.seq[int(pos):]
                self.mseq=seq_m
            else:
                print ("Initial position (%s) not belongs to expected (%s)!")%self.seq[int(pos-1)],self.mutations[pos][0]
                op=raw_input("Do you want to change it anyway? (Y or N)").upper()
                if op=='Y':
                    seq_m+=self.seq[:int(pos)-1]
                    seq_m+=str(self.mutations[pos][1])#aminoacido de substituicao
                    seq_m+=self.seq[int(pos):]
                    self.mseq=str(seq_m)
                elif op=='N':
                    pass 
               
        
        
    def create_mutations(self,pos,am_init,am_fin):
        '''
        Funcao que dadas listas de posicoes, aminoacido inicial e final, 
        as insere na classe para posterior adicao a sequencia normal
        '''
        dic={}
        for p in range(len(pos)):
            dic[pos[p]]=(am_init[p],am_fin[p])#obriga a que as listas sejam dadas no mesmo tamanho
        for key in dic.keys():
            self.mutations[key]=dic[key]
               
        
if __name__=='__main__':

    def test1():
        p=do_seqs()
        #print (p.search_seq_ids("IDH1[Homo Sapiens]"))
        #print p.get_seq("49168486")
        #print p.create_file_with_seq("idh1_normal","49168486")
        #print p.induce_mutation(132,"R","H")
        #print p.create_file_with_seq("IDH1_mutated",0,True)
        #print p.get_gene_file(0,'49456350')
        #print p.convert_into_protein("AATCGA")
        #handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id='49168486')#depois de termos o id, e os valores onde comeca e acaba a seq vamos busca-la
        #read=SeqIO.read(handle,"fasta")
        #print read
        #print p.get_seq(0,'idh1')
        #print p.get_seq('49168486')
        #print p.induce_mutation()
        
        
        #p.create_mutations([2,5,7], ['A','T','C'], ['B','C','T'])
        #print p.search_seq_ids('idh1')
        #print p.get_seq('837785965')
        #print p.get_seq(0, 'idh1')
        #print p.create_file_with_seq('837785965')
        #print p.create_file_with_seq(0,'idh1')
        #print p.convert_into_protein('ATGTGCTGCATGGCA')
        #print p.get_gene_file(0,'3135')
        #print p.get_seq_from_gene('idh1')
    
    
        p.create_file_with_seq('47938312')
        a=raw_input('Posicoes')
        b=raw_input('Inicial')
        a1=[]
        b1=[]
        c1=[]
        c=raw_input('Final')
        for pos in a.split(','):
            a1.append(pos)
        print a1
        for pos in b.split(','):
            b1.append(pos)
        print b1
        for pos in c.split(','):
            c1.append(pos)
        print c1
        p.create_mutations(a1,b1,c1)
        p.induce_mutation()
        p.create_file_with_seq('47938312',0,mutated=True)
        
    
    test1()
    
 
           
        
#handle = Entrez.efetch(db="protein", id="49168486", rettype="gb", retmode="xml")
#results=Entrez.read(handle)
#print results
#Escrever no ficheiro
#SeqIO.write((SeqIO.read(handle,"gb"), 'sequence.gb', "genbank"))

#handle.close() 

#Imprimir na consola os dados obtidos
#record = SeqIO.read("sequence.gb", "genbank")
#print record