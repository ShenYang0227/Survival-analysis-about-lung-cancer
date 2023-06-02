

import json
import os
import shutil

## 1.merge all tsv files in a new folder
def merge_files(file_dir):
    os.mkdir("./gdc")
    for file in os.listdir(file_dir):
        if os.path.isdir(file):
            path = f"{file_dir}/{file}"
            for tsv in os.listdir(path):
                if os.path.splitext(tsv)[-1] == ".tsv":
                    source_path = f"{file_dir}/{file}/{tsv}"
                    target_path = "./gdc"
                    shutil.move(source_path, target_path)

##2.The transformer bwtween file names and TCGA_ID
def Trans_names(fileName):
    trans_name={}
    with open(fileName) as f:
        datas=json.load(f)
        for data in datas:
            file_name=data["file_name"]
            tcga_id = data["associated_entities"][0]["entity_submitter_id"]
            trans_name[file_name]=tcga_id
    return trans_name





##3.count the number of tumor and normal sample
def count_tumor_normal(trans_name):
    normal_list=[]
    tumor_list=[]
    for value in trans_name.values():
        if value.split("-")[3]=="11A":
            normal_list.append(value)
        else:
            tumor_list.append(value)
    num1=len(tumor_list)
    num2=len(normal_list)
    sample_list=tumor_list+normal_list
    return num1,num2,sample_list

##4. obtain the gene matrix and lncRNA matrix

def get_matrix(dir,trans_names,m6A_gene):
    gene_matrix = {}
    lncRNA_matrix={}
    m6A_matrix={}
    id_list=[]
    m6A_list=[]
    with open(m6A_gene) as f:
        for line in f:
            line=line[:-1]
            m6a=line.split("\t")[0]
            m6A_list.append(m6a)

    for tsv in os.listdir(dir):
        for k in trans_names:
            if tsv==k:
                id_list.append(trans_names[k])
                path = f"{dir}/{tsv}"
                with open(path) as f1:
                    for line in f1:
                        line = line[:-1]
                        if line.startswith("ENSG"):
                            data=line.split("\t")
                            gene = data[1]
                            gene_id=data[2]
                            exp = data[-1]
                            if gene not in gene_matrix:
                                gene_matrix.setdefault(gene, []).append(exp)
                            else:
                                gene_matrix[gene].append(exp)
                            if gene_id=="lncRNA":
                                if gene not in lncRNA_matrix:
                                    lncRNA_matrix.setdefault(gene, []).append(exp)
                                else:
                                    lncRNA_matrix[gene].append(exp)
                            for m6A in m6A_list:
                                if m6A==gene:
                                    if gene not in m6A_matrix:
                                        m6A_matrix.setdefault(gene,[]).append(exp)
                                    else:
                                        m6A_matrix.setdefault(gene, []).append(exp)

    return gene_matrix,lncRNA_matrix,m6A_matrix,id_list

##5. write files
def write_files(matrix,file_name,id_list):
    with open(file_name,"w") as fout:
        fout.write("id"+"\t"+"\t".join(id_list)+"\n")
        for k,v in matrix.items():
            fout.write(k + "\t" + "\t".join(v) + "\n")

##6. get survival information
def get_survival(fileName):
    dead_id = []
    alive_id = []
    dead_dict = {}
    alive_dict = {}

    with open(fileName) as f:
        datas = json.load(f)
        for data in datas:
            for k in data.keys():
                if k == "exposures":
                    id = data[k][0]["submitter_id"]
                    id = id.split("_")[0]

                elif k == "demographic":
                    for m in data[k].keys():
                        if m == "vital_status":
                            fustat = data[k][m]
                            if fustat == "Dead":
                                dead_id.append(id)
                                dead_dict[id] = []
                            else:
                                alive_id.append(id)
                                alive_dict[id] = []

                        if m == "days_to_death":
                            time1 = str(data[k][m])
                            dead_dict[id].append(time1)
                            dead_dict[id].append('Dead')

        for data in datas:
            for k in data.keys():
                if k == "exposures":
                    id = data[k][0]["submitter_id"]
                    id = id.split("_")[0]

                    for i in alive_id:
                        if id == i:
                            for k in data.keys():
                                if k == "diagnoses":
                                    time2 = str(data[k][0]["days_to_last_follow_up"])

                                    alive_dict[id].append(time2)
                                    alive_dict[id].append("Alive")

        surv_infor = dict(list(alive_dict.items()) + list(dead_dict.items()))

        for k, v in surv_infor.items():
            if v == []:
                surv_infor[k].append("NA")
                surv_infor[k].append("NA")

        with open("time.txt", "w") as fin:
            fin.write("id" + "\t" + "futime" + "\t" + "fustat" + "\n")
            for k, v in surv_infor.items():
                fin.write(k+"\t"+"\t".join(v)+"\n")

##7. get the clinical information
def get_clinical(file_name):
    with open(file_name) as f1:
        clin_infor = {}
        datas = json.load(f1)
        for data in datas:
            for k in data.keys():
                if k == "exposures":
                    id = data[k][0]["submitter_id"]
                    id = id.split("_")[0]
                    clin_infor[id] = []

                if k == "demographic":
                    gender = data[k]["gender"]
                    age = data[k]["age_at_index"]
                    clin_infor[id].append(gender)
                    clin_infor[id].append(str(age))

                if k == "diagnoses":
                    for m, n in data[k][0].items():
                        if m == "ajcc_pathologic_stage":
                            stage = data[k][0][m]
                            clin_infor[id].append(stage)
                        if m == "ajcc_pathologic_t":
                            T = data[k][0][m]
                            clin_infor[id].append(T)

                        if m == "ajcc_pathologic_n":
                            N = data[k][0][m]
                            clin_infor[id].append(N)
                        if m == "ajcc_pathologic_m":
                            M = data[k][0][m]
                            clin_infor[id].append(M)

        with open("clinical.txt", "w") as fout:
            fout.write("id" + "\t"+"stage" + "\t" + "T" + "\t" + "N" + "\t" + "M" + "\t" + "gender" + "\t" + "age" + "\n")
            for k, v in clin_infor.items():
                fout.write(k + "\t" + "\t".join(v) + "\n")





if __name__=="__main__":
    trans=Trans_names("./metadata.cart.2023-05-10.json")
    tumor_num,normal_num,sample_id=count_tumor_normal(trans)
    print("the number of tumor:",tumor_num)
    print("the number of normal:",normal_num)
    gene_matrix,lncRNA_matrix,m6A_matrix,id_list=get_matrix("./gdc",trans,"./gene.txt")
    write_files(gene_matrix,"./gene_matrix.txt",id_list)
    write_files(lncRNA_matrix, "./lncRNA_matrix.txt", id_list)
    write_files(m6A_matrix, "./m6A_matrix.txt", id_list)
    get_survival("./clinical.cart.2023-05-10.json")
    get_clinical("./clinical.cart.2023-05-10.json")













