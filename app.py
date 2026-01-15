from flask import Flask, render_template, request
from Bio import Entrez
from datetime import datetime, timedelta
import requests
import os

# Configuração do Flask
base_dir = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, template_folder=os.path.join(base_dir, 'templates'))

# --- CONFIGURAÇÕES ---
Entrez.email = "seu_email@email.com"
# A sua chave de API
GEMINI_API_KEY = "AIzaSyASnfSyvIrKmPKj2VHt4YOY3Vcfh6Vs_g0"

def gerar_conteudo_gemini(prompt):
    """Função robusta para chamar a API do Gemini utilizando a versão v1 estável."""
    # Alterado para v1 para evitar problemas de compatibilidade com modelos flash
    url = f"https://generativelanguage.googleapis.com/v1/models/gemini-1.5-flash:generateContent?key={GEMINI_API_KEY}"
    
    payload = {
        "contents": [{
            "parts": [{"text": prompt}]
        }]
    }
    
    headers = {
        "Content-Type": "application/json"
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers, timeout=30)
        data = response.json()
        
        # Verificação da resposta
        if 'candidates' in data and len(data['candidates']) > 0:
            return data['candidates'][0]['content']['parts'][0]['text']
        elif 'error' in data:
            return f"Erro da API Gemini: {data['error'].get('message', 'Erro desconhecido')}"
        else:
            return "A IA não conseguiu gerar uma resposta. Tente novamente."
            
    except Exception as e:
        return f"Erro de conexão com a IA: {str(e)}"

def processar_termos_busca(termos_usuario):
    """Traduz termos para inglês técnico."""
    prompt = f"Traduza para termos científicos MeSH em inglês, retorne APENAS os termos entre parênteses: {termos_usuario}"
    resultado = gerar_conteudo_gemini(prompt)
    
    # Limpeza básica caso a IA responda com texto extra
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
        
        # 1. Tradução
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
                
                # 3. Resumo Clínico focado em Reabilitação e Hipertrofia
                prompt_resumo = f"""
                Analise estes artigos científicos sobre {tema_original}.
                Traduza e resuma para Português com foco em:
                1. Conclusões principais.
                2. Aplicação prática para Recuperação de Lesões.
                3. Aplicação prática para Hipertrofia.
                4. Links de referência: https://pubmed.ncbi.nlm.nih.gov/ID/
                
                ARTIGOS:
                {texto_abstracts}
                """
                relatorio = gerar_conteudo_gemini(prompt_resumo)
            else:
                relatorio = "Nenhum artigo recente encontrado. Tente termos mais amplos ou um período maior."
        except Exception as e:
            relatorio = f"Erro no sistema: {str(e)}"
            
    return render_template('index.html', relatorio=relatorio, tema=tema_exibicao, periodo=periodo_selecionado)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
