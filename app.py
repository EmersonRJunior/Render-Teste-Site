from flask import Flask, render_template, request
from Bio import Entrez
from datetime import datetime, timedelta
import google.generativeai as genai
import os
import re

# Configuração do Flask
base_dir = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, template_folder=os.path.join(base_dir, 'templates'))

# --- CONFIGURAÇÕES ---
# Coloque um email válido para o PubMed não bloquear o seu IP
Entrez.email = "seu_email@exemplo.com"
GEMINI_API_KEY = "AIzaSyASnfSyvIrKmPKj2VHt4YOY3Vcfh6Vs_g0"

# Configuração Oficial da API Google
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel('gemini-1.5-flash')

def gerar_conteudo_gemini(prompt):
    """Gera conteúdo usando a SDK oficial do Google."""
    try:
        response = model.generate_content(prompt)
        if response and response.text:
            return response.text
        return ""
    except Exception as e:
        print(f"Erro Gemini: {e}")
        return ""

def limpar_query(texto):
    """Extrai apenas o que está entre parênteses ou limpa lixo da resposta da IA."""
    # Tenta encontrar conteúdo entre parênteses
    match = re.search(r'\((.*?)\)', texto)
    if match:
        return match.group(0)
    # Se não houver parênteses, remove caracteres especiais e mantém palavras-chave
    limpo = re.sub(r'[^a-zA-Z0-9\s"()]', '', texto)
    return limpo.strip()

def processar_termos_busca(termos_usuario):
    """Traduz termos do utilizador para inglês científico MeSH."""
    prompt = (
        f"Atue como um bibliotecário científico. Traduza o seguinte tema de fisioterapia/treino "
        f"para uma query de busca PubMed (MeSH terms) em inglês. "
        f"Retorne APENAS a query entre parênteses. Exemplo: (Hypertrophy AND Resistance Training). "
        f"Tema: {termos_usuario}"
    )
    resultado = gerar_conteudo_gemini(prompt)
    query = limpar_query(resultado)
    return query if query else f"({termos_usuario})"

def buscar_detalhes_pubmed(id_list):
    """Recupera os abstracts do PubMed com tratamento de erro."""
    if not id_list: return ""
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="text")
        dados = handle.read()
        handle.close()
        return dados
    except Exception as e:
        print(f"Erro ao buscar abstracts: {e}")
        return ""

@app.route('/', methods=['GET', 'POST'])
def index():
    relatorio = None
    tema_exibicao = ""
    query_executada = ""
    
    if request.method == 'POST':
        tema_original = request.form.get('tema')
        periodo = request.form.get('periodo', '30')
        tema_exibicao = tema_original
        
        # 1. Tradução e Limpeza
        query_ingles = processar_termos_busca(tema_original)
        
        # 2. Formatação da Data (PubMed format: YYYY/MM/DD)
        data_limite = (datetime.now() - timedelta(days=int(periodo))).strftime("%Y/%m/%d")
        query_final = f"{query_ingles} AND (\"{data_limite}\"[Date - Publication] : \"3000\"[Date - Publication])"
        query_executada = query_final # Para debug visual se necessário
        
        try:
            # Busca no PubMed
            handle = Entrez.esearch(db="pubmed", term=query_final, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            
            # Se não encontrar nada com data, tenta busca geral (sem data) para não dar vazio
            if not ids:
                handle = Entrez.esearch(db="pubmed", term=query_ingles, retmax=3)
                record = Entrez.read(handle)
                handle.close()
                ids = record.get("IdList", [])

            if ids:
                texto_abstracts = buscar_detalhes_pubmed(ids)
                
                # 3. Resumo com foco em Recuperação e Hipertrofia
                prompt_resumo = f"""
                Analise estes resumos científicos sobre {tema_original}:
                {texto_abstracts}
                
                Escreva um relatório em PORTUGUÊS (Brasil) focado em:
                - Principais achados científicos.
                - Aplicação prática para REABILITAÇÃO/RECUPERAÇÃO.
                - Aplicação prática para HIPERTROFIA.
                - Liste os links: https://pubmed.ncbi.nlm.nih.gov/ID/
                """
                relatorio = gerar_conteudo_gemini(prompt_resumo)
            else:
                relatorio = f"Não encontramos artigos recentes para: {query_ingles}. Tente termos mais simples como 'ACL recovery' ou 'Muscle Hypertrophy'."
        except Exception as e:
            relatorio = f"Erro no processo de busca: {str(e)}"
            
    return render_template('index.html', relatorio=relatorio, tema=tema_exibicao)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)