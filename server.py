import math
import os
from flask import Flask, jsonify, request,send_from_directory
from flask_cors import CORS
from pymongo import MongoClient

app = Flask(__name__)

# 启用 CORS
CORS(app)

# 配置 MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client['gene_data']  # 替换为您的数据库名

def replace_nan_with_dash(record):
    for key, value in record.items():
        if isinstance(value, float) and math.isnan(value):
            record[key] = '-'
        elif isinstance(value, dict):
            replace_nan_with_dash(value)
    return record


@app.route('/search', methods=['GET'])
def search():
    term = request.args.get('query')
    if term:
        # 在 'zju' 集合中搜索 GeneID 包含 term 的记录
        records = db.zju.find({"GeneID": {"$regex": term}}, {"GeneID": 1}).limit(5)
        # 获取不重复的 GeneID
        gene_ids = {record['GeneID'] for record in records}
        return jsonify(list(gene_ids))
    else:
        return jsonify([])
@app.route('/search2', methods=['GET'])
def search2():
    term = request.args.get('query')
    if term:
        # 在 'zju' 集合中搜索 GeneID 包含 term 的记录
        records = db.zju.find({"Chr ID": {"$regex": term}}, {"Chr ID": 1}).limit(5)
        # 获取不重复的 GeneID
        gene_ids = {record['Chr ID'] for record in records}
        return jsonify(list(gene_ids))
    else:
        return jsonify([])

@app.route('/get-gene-info', methods=['GET'])
def get_gene_info():
    gene_id = request.args.get('geneId')
    if gene_id:
        record = db.zju.find_one({"GeneID": gene_id}, {"_id": 0})
        if record:
            record = replace_nan_with_dash(record)
        return jsonify(record) if record else jsonify({})
    else:
        return jsonify({})

@app.route('/get-gene-info2', methods=['GET'])
def get_gene_info2():
    gene_id = request.args.get('chrid')
    if gene_id:
        record = db.zju.find_one({"Chr ID": gene_id}, {"_id": 0})
        if record:
            record = replace_nan_with_dash(record)
        return jsonify(record) if record else jsonify({})
    else:
        return jsonify({})

@app.route('/chartdata/<geneId>', methods=['GET'])
def get_chart_data(geneId):
    geneId = 'MrScaffold_086G1'
    records = db.sample6_exp.find({'Gene ID': geneId}, {'_id': 0})
    data = list(records)

    # 返回JSON响应
    return jsonify(data)


@app.route('/searchvarian', methods=['POST'])
def searchvarian():
    data = request.json
    # db.snp.create_index([("ID", 1)])
    search_text = data.get('searchText3')
    results = db.snp.find({"ID": {"$regex": search_text}},{"_id":0})
    results_list = [doc for doc in results]
    return jsonify(results_list)


@app.route('/getvarianoption', methods=['GET'])
def get_varian_option():
    distinct_chroms = db.snp.distinct('#CHROM')
    chrom_stats = {}

    for chrom in distinct_chroms:
        # 聚合查询以获取每个chrom的最大和最小POS值
        aggregation = [
            {"$match": {"#CHROM": chrom}},
            {"$group": {
                "_id": "$#CHROM",
                "minPOS": {"$min": "$POS"},
                "maxPOS": {"$max": "$POS"}
            }}
        ]
        result = list(db.snp.aggregate(aggregation))
        if result:
            chrom_stats[str(chrom)] = {
                "min": result[0]['minPOS'],
                "max": result[0]['maxPOS']
            }

    return jsonify({"chrom": distinct_chroms, **chrom_stats})

@app.route('/getvairanFields', methods=['GET'])
def getvairanFields():
    res = db.snp.find_one()
    m_fields = [key for key in res if key.startswith('M.')]
    return jsonify(m_fields)

@app.route('/gettranscriptomefields', methods=['GET'])
def gettranscriptomefields():
    blist=['_id','Gene ID']
    res = db.sample6_exp.find_one()
    m_fields = [key for key in res if key not in blist]
    return jsonify(m_fields)

@app.route('/multivarian', methods=['POST'])
def query_snp():
    data = request.json
    selected_option = int(data['selectedOption'])
    varian_start = int(data['varianstart'])
    varian_end = int(data['varianend'])
    rightitems = set(data['rightItems'])  # 使用集合确保唯一性和快速查找

    collection = db.snp
    query_result = collection.find({
        "#CHROM": selected_option,
        "POS": {"$gte": varian_start, "$lte": varian_end}
    }, {'_id': 0})

    # 处理查询结果
    result_list = []
    for record in query_result:
        # 移除不在 rightitems 中的 M 开头的字段
        filtered_record = {k: v for k, v in record.items() if not (k.startswith('M.') and k not in rightitems)}

        # 计算符合条件的 M 开头字段数量
        count = sum(1 for k in rightitems if k in filtered_record and '1' in filtered_record[k])

        # 将 count 值添加到记录中
        filtered_record['count'] = str(count)+'/'+str(len(rightitems) if rightitems else 0)
        result_list.append(filtered_record)

    return jsonify(result_list)

@app.route('/transcriptome_data', methods=['POST'])
def transcriptome_data():
    data = request.json
    idlist = data.get('idlist', [])
    rightitems = set(data.get('rightitems', []))

    query_conditions = {'Gene ID': {'$in': idlist}}

    # 根据 rightitems 的内容构建投影
    if rightitems:
        # rightitems 非空，只包含指定字段
        projection = {item: 1 for item in rightitems}
        projection['_id'] = 0
        projection['Gene ID'] = 1
    else:
        # rightitems 为空，包含所有字段除了 '_id'
        projection = {'_id': 0}

    record = db.sample6_exp.find(query_conditions, projection)
    result_list = [replace_nan_with_dash(rec) for rec in record]

    return jsonify(result_list)


@app.route('/getmetabolomicsfields', methods=['GET'])
def getmetabolomicsfields():
    blist=['_id',"Ion_mode","Q1","Q3","RT","Compound_name","Class_I","Class_II",'ID']
    res = db.meta.find_one()
    m_fields = [key for key in res if key not in blist]
    return jsonify(m_fields)

@app.route('/metabolomics_data', methods=['POST'])
def metabolomics_data():
    data = request.json
    idlist = data.get('idlist', [])
    rightitems = set(data.get('rightitems', []))

    # 构建查询条件
    query_conditions = {'Compound_name': {'$in': idlist}}

    # 从 MongoDB 获取数据（只排除 '_id' 字段）
    records = db.meta.find(query_conditions, {'_id': 0})

    # 默认需要保留的字段集合
    default_fields = {'ID', 'Ion_mode', 'Q1', 'Q3', 'RT', 'Compound_name', 'Class_I', 'Class_II'}

    # 处理结果
    result_list = []
    for record in records:
        # 只保留 rightitems 中的字段和默认字段
        filtered_record = {key: value for key, value in record.items() if key in rightitems or key in default_fields}
        result_list.append(replace_nan_with_dash(filtered_record))

    return jsonify(result_list)

files_directory = '../files'
@app.route('/get-files')
def get_files():
    files = os.listdir(files_directory)
    data = [{'name': f, 'url': f'/download/{f}'} for f in files]
    return jsonify(data)

@app.route('/download/<filename>')
def download_file(filename):
    print(filename)
    return send_from_directory(files_directory, filename, as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
