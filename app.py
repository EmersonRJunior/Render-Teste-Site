from flask import Flask, render_template, request, flash
from Bio import Entrez
import google.generativeai as genai
import os
from dotenv import load_dotenv

load_dotenv()

app = Flask(__name__)
app.secret_key = os.getenv("FLASK_SECRET_KEY", "chave-secreta-fisioterapia")

# Configurações
Entrez.email = os.getenv("PUBMED_EMAIL", "seu-email@exemplo.com")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")

if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)
    model = genai.GenerativeModel('gemini-1.5-flash')

@app.route('/', methods=['GET', 'POST'])
def index():
    resultados = None
    tema_exibicao = ""
    
    if request.method == 'POST':
        tema = request.form.get('tema', '').strip()
        tema_exibicao = tema
        
        try:
            # 1. Busca direta no PubMed usando o termo do usuário
            # Buscamos os 5 IDs mais relevantes
            with Entrez.esearch(db="pubmed", term=tema, retmax=5) as handle:
                record = Entrez.read(handle)
                ids = record.get("IdList", [])

            if ids:
                # 2. Busca os títulos dos artigos encontrados
                with Entrez.efetch(db="pubmed", id=ids, rettype="pref", retmode="xml") as handle:
                    detalhes = Entrez.read(handle)
                
                lista_artigos = []
                for artigo in detalhes['PubmedArticle']:
                    titulo = artigo['MedlineCitation']['Article']['ArticleTitle']
                    pmid = artigo['MedlineCitation']['PMID']
                    link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    lista_artigos.append(f"Título: {titulo}\nLink: {link}")

                # 3. IA apenas organiza a lista de forma amigável
                texto_para_ia = "\n\n".join(lista_artigos)
                prompt = f"Abaixo estão artigos sobre '{tema}'. Formate-os em uma lista organizada em português: \n\n{texto_para_ia}"
                
                response = model.generate_content(prompt)
                resultados = response.text
            else:
                resultados = "Nenhum artigo encontrado para este tema específico."

        except Exception as e:
            print(f"Erro: {e}")
            flash("Erro ao conectar com as bases de dados.")
            
    return render_template('index.html', relatorio=resultados, tema=tema_exibicao)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
