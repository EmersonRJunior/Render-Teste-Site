from flask import Flask, render_template, request
from Bio import Entrez
from datetime import datetime, timedelta
import google.generativeai as genai
import os

# Configuração do Flask
base_dir = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, template_folder=os.path.join(base_dir, 'templates'))

# --- CONFIGURAÇÕES ---
Entrez.email = "seu_email@email.com"
GEMINI_API_KEY = "AIzaSyASnfSyvIrKmPKj2VHt4YOY3Vcfh6Vs_g0"

# Configuração Oficial da API Google
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel('gemini-1.5-flash')

def gerar_conteudo_gemini(prompt):
    """Usa a biblioteca oficial para gerar conteúdo, mais estável no Render."""
    try:
        response = model.generate_content(prompt)
        if response and response.text:
            return response.text
        return "A IA não retornou dados. Tente novamente."
    except Exception as e:
        return f"Erro na API Gemini (Região/Chave): {str(e)}"

def processar_termos_busca(termos_usuario):
    """Traduz termos para inglês técnico MeSH."""
    prompt = f"Traduza os termos de reabilitação/fisioterapia para MeSH terms em inglês. Retorne APENAS os termos entre parênteses: {termos_usuario}"
    resultado = gerar_conteudo_gemini(prompt)
    
    if "(" in resultado and ")" in resultado:
        return resultado[resultado.find("("):resultado.find(")")+1]
    return f"({termos_usuario})"

def buscar_detalhes_pubmed(id_list):
    """Recupera os abstracts do PubMed."""
    if not id_list: return ""
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="text")
        dados = handle.read()
        handle.close()
        return dados
    except:
        return ""

@app.route('/', methods=['GET', 'POST'])
def index():
    relatorio = None
    tema_exibicao = ""
    periodo_selecionado = "30"
    
    if request.method == 'POST':
        tema_original = request.form.get('tema')
        periodo_selecionado = request.form.get('periodo', '30')
        tema_exibicao = tema_original
        
        # 1. Tradução via IA
        query_ingles = processar_termos_busca(tema_original)
        
        # 2. Busca PubMed
        data_inicio = (datetime.now() - timedelta(days=int(periodo_selecionado))).strftime("%Y/%m/%d")
        query_final = f"{query_ingles} AND ({data_inicio}[PDAT] : 3000[PDAT])"
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query_final, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            
            if ids:
                texto_abstracts = buscar_detalhes_pubmed(ids)
                
                # 3. Resumo com foco em Recuperação e Hipertrofia
                prompt_resumo = f"""
                Você é um especialista em fisiologia do exercício. Analise estes abstracts sobre {tema_original}:
                {texto_abstracts}
                
                Gere um relatório em PORTUGUÊS com:
                - Resumo das evidências atuais.
                - Aplicação prática para Recuperação de Lesões.
                - Aplicação prática para Hipertrofia.
                - Links (https://pubmed.ncbi.nlm.nih.gov/ID/).
                """
                relatorio = gerar_conteudo_gemini(prompt_resumo)
            else:
                relatorio = "Nenhum artigo encontrado para este período."
        except Exception as e:
            relatorio = f"Erro no servidor: {str(e)}"
            
    return render_template('index.html', relatorio=relatorio, tema=tema_exibicao, periodo=periodo_selecionado)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
