from flask import Flask, render_template, request
from Bio import Entrez
from datetime import datetime, timedelta
import requests
import json
import os

# Configuração para encontrar a pasta templates no Render
base_dir = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, template_folder=os.path.join(base_dir, 'templates'))

# --- CONFIGURAÇÕES ---
# Identificação para o PubMed
Entrez.email = "seu_email@email.com"

# Tua Chave de API do Gemini
GEMINI_API_KEY = "AIzaSyASnfSyvIrKmPKj2VHt4YOY3Vcfh6Vs_g0"

def obter_modelo_disponivel():
    """Consulta a API do Google para encontrar um modelo ativo (v1beta)."""
    url = f"https://generativelanguage.googleapis.com/v1beta/models?key={GEMINI_API_KEY}"
    try:
        response = requests.get(url, timeout=10)
        data = response.json()
        if 'models' in data:
            for m in data['models']:
                if 'generateContent' in m.get('supportedGenerationMethods', []):
                    return m['name']
        return "models/gemini-1.5-flash"
    except Exception:
        return "models/gemini-1.5-flash"

def processar_termos_busca(termos_usuario):
    """Traduz termos do português para inglês técnico e formata a query do PubMed."""
    modelo_path = obter_modelo_disponivel()
    url = f"https://generativelanguage.googleapis.com/v1beta/{modelo_path}:generateContent?key={GEMINI_API_KEY}"
    
    prompt = f"Traduza os seguintes termos de fisioterapia/reabilitação para termos científicos em inglês (MeSH). Retorne apenas os termos em inglês separados por ' AND ' e entre parênteses. Entrada: {termos_usuario}"
    
    payload = {"contents": [{"parts": [{"text": prompt}]}]}
    try:
        response = requests.post(url, json=payload, timeout=20)
        result = response.json()
        query = result['candidates'][0]['content']['parts'][0]['text'].replace('`', '').strip()
        return query
    except Exception:
        # Se falhar a tradução, faz um replace básico para não travar a busca
        return f"({termos_usuario.replace(',', ' AND ')})"

def buscar_detalhes_pubmed(id_list):
    """Recupera os abstracts do PubMed."""
    if not id_list: return ""
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="text")
        dados = handle.read()
        handle.close()
        return dados
    except Exception:
        return ""

def gerar_resumo_ia(texto_artigos, tema_original):
    """Gera o relatório clínico final em Português."""
    if not texto_artigos: return "Dados insuficientes para análise."
    
    modelo_path = obter_modelo_disponivel()
    url = f"https://generativelanguage.googleapis.com/v1beta/{modelo_path}:generateContent?key={GEMINI_API_KEY}"
    
    prompt = f"""
    Você é um especialista em Fisiologia e Reabilitação com foco em Recuperação de Lesões e Hipertrofia.
    Analise os estudos sobre: {tema_original}.
    Traduza as conclusões para PORTUGUÊS técnico.

    ESTRUTURA DO RELATÓRIO:
    1. RESUMO DAS EVIDÊNCIAS: O que os estudos mais recentes dizem?
    2. RECUPERAÇÃO DE LESÕES: Implicações práticas para reabilitação.
    3. HIPERTROFIA: Implicações para o ganho de massa muscular.
    4. PROTOCOLO SUGERIDO: Aplicação prática imediata.

    Inclua os links no final: https://pubmed.ncbi.nlm.nih.gov/ID/
    
    TEXTO DOS ARTIGOS:
    {texto_artigos}
    """
    
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "safetySettings": [
            {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
        ]
    }
    
    try:
        response = requests.post(url, json=payload, timeout=45)
        result = response.json()
        return result['candidates'][0]['content']['parts'][0]['text']
    except Exception as e:
        return f"Erro ao gerar resumo: {str(e)}"

@app.route('/', methods=['GET', 'POST'])
def index():
    relatorio = None
    tema_exibicao = ""
    periodo = "30"
    
    if request.method == 'POST':
        tema_original = request.form.get('tema')
        periodo = request.form.get('periodo', '30')
        tema_exibicao = tema_original
        
        # 1. Tradução e Expansão via IA
        query_ingles = processar_termos_busca(tema_original)
        
        # 2. Busca Temporal
        data_inicio = (datetime.now() - timedelta(days=int(periodo))).strftime("%Y/%m/%d")
        query_final = f"{query_ingles} AND ({data_inicio}[PDAT] : 3000[PDAT])"
        
        try:
            handle = Entrez.esearch(db="pubmed", term=query_final, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            
            if ids:
                texto_abstracts = buscar_detalhes_pubmed(ids)
                relatorio = gerar_resumo_ia(texto_abstracts, tema_original)
            else:
                relatorio = f"Nenhum artigo encontrado para '{tema_original}' nos últimos {periodo} dias."
        except Exception as e:
            relatorio = f"Erro no servidor: {str(e)}"
            
    return render_template('index.html', relatorio=relatorio, tema=tema_exibicao, periodo=periodo)

if __name__ == '__main__':
    # Configuração vital para o Render (usa a porta definida pela variável de ambiente)
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)