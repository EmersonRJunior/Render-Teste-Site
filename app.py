from flask import Flask, render_template, request
from Bio import Entrez
from datetime import datetime, timedelta
import google.generativeai as genai
import os
import re

app = Flask(__name__)

# --- CONFIGURAÇÕES ---
Entrez.email = "seu_email@exemplo.com"
# Recomenda-se setar no terminal: export GEMINI_API_KEY='sua_chave'
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY", "SUA_CHAVE_AQUI_CASO_NAO_USE_ENV")

genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel('gemini-1.5-flash')

def gerar_conteudo_gemini(prompt):
    try:
        response = model.generate_content(prompt)
        return response.text if response else ""
    except Exception as e:
        print(f"Erro Gemini: {e}")
        return ""

def limpar_query(texto):
    # Procura conteúdo dentro de parênteses primeiro
    match = re.search(r'\((.*?)\)', texto)
    if match:
        return f"({match.group(1)})"
    
    # Se não achar, remove apenas aspas e limpa espaços extras
    limpo = texto.replace('"', '').replace('`', '').strip()
    return limpo

def processar_termos_busca(termos_usuario):
    prompt = (
        f"Converta o tema '{termos_usuario}' em uma query de busca PubMed profissional (MeSH terms). "
        "Retorne APENAS a query entre parênteses. Exemplo: (Resistance Training AND Muscle Hypertrophy)."
    )
    resultado = gerar_conteudo_gemini(prompt)
    return limpar_query(resultado)

def buscar_detalhes_pubmed(id_list):
    if not id_list: return ""
    try:
        # efetch funciona melhor com IDs separados por vírgula
        ids_string = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", id=ids_string, rettype="abstract", retmode="text")
        dados = handle.read()
        handle.close()
        return dados
    except Exception as e:
        print(f"Erro ao buscar detalhes: {e}")
        return ""

@app.route('/', methods=['GET', 'POST'])
def index():
    relatorio = None
    tema_exibicao = ""
    
    if request.method == 'POST':
        tema_original = request.form.get('tema')
        periodo = request.form.get('periodo', '30')
        tema_exibicao = tema_original
        
        # 1. Tradução
        query_ingles = processar_termos_busca(tema_original)
        
        # 2. Busca no PubMed
        try:
            data_limite = (datetime.now() - timedelta(days=int(periodo))).strftime("%Y/%m/%d")
            query_final = f"{query_ingles} AND (\"{data_limite}\"[Date - Publication] : \"3000\"[Date - Publication])"
            
            handle = Entrez.esearch(db="pubmed", term=query_final, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])

            # Fallback se não houver artigos recentes
            if not ids:
                handle = Entrez.esearch(db="pubmed", term=query_ingles, retmax=3)
                record = Entrez.read(handle)
                handle.close()
                ids = record.get("IdList", [])

            if ids:
                texto_abstracts = buscar_detalhes_pubmed(ids)
                
                # Criar a lista de links reais para passar ao Gemini
                links_reais = "\n".join([f"https://pubmed.ncbi.nlm.nih.gov/{id_}/" for id_ in ids])
                
                prompt_resumo = f"""
                Analise estes abstracts:
                {texto_abstracts}
                
                Escreva um relatório em PORTUGUÊS para o tema '{tema_original}':
                - Resumo dos achados.
                - Aplicação em REABILITAÇÃO.
                - Aplicação em HIPERTROFIA.
                
                Use estes links reais para as referências:
                {links_reais}
                """
                relatorio = gerar_conteudo_gemini(prompt_resumo)
            else:
                relatorio = "Nenhum artigo encontrado para este tema."
                
        except Exception as e:
            relatorio = f"Erro técnico: {str(e)}"
            
    return render_template('index.html', relatorio=relatorio, tema=tema_exibicao)

if __name__ == '__main__':
    app.run(debug=True) # debug=True ajuda a ver erros no navegador
