from flask import Flask, render_template, request, flash
from Bio import Entrez
from datetime import datetime, timedelta
import google.generativeai as genai
import os
import re
from dotenv import load_dotenv

# Carrega variáveis do arquivo .env
load_dotenv()

app = Flask(__name__)
app.secret_key = os.getenv("FLASK_SECRET_KEY", "chave-secreta-padrao-para-dev")

# --- CONFIGURAÇÕES VIA AMBIENTE ---
# Nunca deixe sua API Key exposta no código ao subir para o GitHub
Entrez.email = os.getenv("PUBMED_EMAIL", "seu-email@exemplo.com")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")

if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)
    model = genai.GenerativeModel('gemini-1.5-flash')
else:
    print("AVISO: GEMINI_API_KEY não encontrada. O sistema de IA não funcionará.")

def gerar_conteudo_gemini(prompt):
    """Gera conteúdo com tratamento de erro e segurança."""
    if not GEMINI_API_KEY:
        return "Erro: API Key não configurada."
    try:
        response = model.generate_content(prompt)
        return response.text if response and response.text else ""
    except Exception as e:
        app.logger.error(f"Erro na API Gemini: {e}")
        return ""

def processar_termos_busca(termos_usuario):
    """Converte termos do usuário para MeSH utilizando IA."""
    # Limpeza básica contra injeção de prompt
    termos_limpos = re.sub(r'[^\w\s]', '', termos_usuario)
    
    prompt = (
        f"Atue como um bibliotecário científico. Traduza o tema abaixo para uma "
        f"query PubMed (MeSH terms) em inglês. Retorne APENAS a query entre parênteses. "
        f"Tema: {termos_limpos}"
    )
    
    resultado = gerar_conteudo_gemini(prompt)
    
    # Extração robusta da query entre parênteses
    match = re.search(r'\((.*?)\)', resultado)
    return match.group(0) if match else f"({termos_limpos})"

def buscar_detalhes_pubmed(id_list):
    """Busca abstracts de forma segura."""
    if not id_list: return ""
    try:
        with Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="text") as handle:
            return handle.read()
    except Exception as e:
        app.logger.error(f"Erro PubMed Fetch: {e}")
        return ""

@app.route('/', methods=['GET', 'POST'])
def index():
    relatorio = None
    tema_exibicao = ""
    
    if request.method == 'POST':
        tema_original = request.form.get('tema', '').strip()
        periodo = request.form.get('periodo', '30')
        
        if not tema_original:
            flash("Por favor, digite um tema para pesquisa.")
            return render_template('index.html')

        tema_exibicao = tema_original
        query_ingles = processar_termos_busca(tema_original)
        
        # Formatação de data PubMed
        data_limite = (datetime.now() - timedelta(days=int(periodo))).strftime("%Y/%m/%d")
        query_final = f"{query_ingles} AND (\"{data_limite}\"[Date - Publication] : \"3000\"[Date - Publication])"
        
        try:
            # Busca IDs
            with Entrez.esearch(db="pubmed", term=query_final, retmax=5) as handle:
                record = Entrez.read(handle)
                ids = record.get("IdList", [])

            # Fallback caso não encontre artigos recentes
            if not ids:
                with Entrez.esearch(db="pubmed", term=query_ingles, retmax=3) as handle:
                    record = Entrez.read(handle)
                    ids = record.get("IdList", [])

            if ids:
                texto_abstracts = buscar_detalhes_pubmed(ids)
                prompt_resumo = (
                    f"Analise estes resumos científicos sobre {tema_original}:\n{texto_abstracts}\n\n"
                    f"Escreva um relatório em PORTUGUÊS (Brasil) focado em:\n"
                    f"- Principais achados.\n- Aplicação em REABILITAÇÃO.\n- Aplicação em HIPERTROFIA.\n"
                    f"Ao final, liste os links no formato: https://pubmed.ncbi.nlm.nih.gov/ID/"
                )
                relatorio = gerar_conteudo_gemini(prompt_resumo)
            else:
                relatorio = "Nenhum artigo encontrado para este tema."
                
        except Exception as e:
            app.logger.error(f"Erro Geral: {e}")
            relatorio = "Ocorreu um erro técnico ao processar sua busca."
            
    return render_template('index.html', relatorio=relatorio, tema=tema_exibicao)

if __name__ == '__main__':
    # Configuração para deploy (Heroku, Render, etc)
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
