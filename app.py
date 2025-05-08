from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

app = Flask(__name__)

@app.route('/generate-name', methods=['POST'])
def generate_name():
    data = request.get_json()
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({'error': 'No SMILES provided'}), 400

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return jsonify({'error': 'Invalid SMILES'}), 400

        # 예시: 분자식 계산 (IUPAC 이름 생성은 추후 Open Babel 또는 PubChem으로)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        return jsonify({'formula': formula})

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)
