"""Format retained simulation rows as complete tables grouped by sample size.
No results are recomputed. Run after the original simulation-table generator.
"""
from pathlib import Path
import argparse,re,json

def braced(text,start):
    assert text[start]=='{';depth=0
    for j in range(start,len(text)):
        if text[j]=='{':depth+=1
        elif text[j]=='}':
            depth-=1
            if depth==0:return text[start+1:j],j+1
    raise ValueError('Unbalanced braces')

def split_tables(text):
    audit=[]
    def convert(match):
        block=match.group(0);spec,end=braced(block,block.index('{',len(r'\begin{longtable}')))
        capstart=block.index(r'\caption')+len(r'\caption');caption,_=braced(block,capstart)
        label=re.search(r'\\label\{([^}]+)\}',block).group(1)
        header=block[block.index(r'\toprule'):block.index(r'\endfirsthead')].strip()
        body=block.split(r'\endlastfoot',1)[1].split(r'\end{longtable}',1)[0]
        rows=[line.strip() for line in body.splitlines() if re.match(r'^\s*(20|33|60|120)\s*&',line)]
        assert len(rows)==108,(label,len(rows))
        parts=[];recovered=[]
        for n in [20,33,60,120]:
            subset=[line for line in rows if re.match(fr'^{n}\s*&',line)];assert len(subset)==27
            recovered.extend(subset)
            caption_n=caption.replace('Complete 108-cell repeated-summer results:',f'Simulation results for $R={n}$ summers:')
            label_n=label if n==20 else label+f'-r{n}'
            parts.append(r'''\begin{table}[H]
\centering
\scriptsize
\caption{'''+caption_n+'}\n'+r'\label{'+label_n+'}\n'+r'\begin{tabular}{'+spec+'}\n'+header+'\n'+'\n'.join(subset)+'\n'+r'''\bottomrule
\end{tabular}
\end{table}
''')
        assert rows==recovered
        audit.append(dict(label=label,original_rows=108,parts=4,rows_per_part=27,rows_preserved_in_order=True))
        return '\n'.join(parts)
    result=re.sub(r'\\begin\{longtable\}.*?\\end\{longtable\}',convert,text,flags=re.S)
    if audit:
        # Keep the small table type local to this generated input.
        result='% Simulation tables grouped by R; all retained rows preserved.\n\\begingroup\n'+result.rstrip()+'\n\\endgroup\n'
    return result,audit

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);parser.add_argument('--audit',type=Path)
    args=parser.parse_args();source=args.input.read_text();text,audit=split_tables(source)
    args.output.write_text(text)
    if args.audit:args.audit.write_text(json.dumps(audit,indent=2))
    print('Formatted',len(audit),'long tables; retained',sum(x['original_rows'] for x in audit),'rows.')
