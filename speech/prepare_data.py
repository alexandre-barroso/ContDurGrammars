"""Parse public ARCTIC Festival relations; no acoustic feature enters a predictor."""
from pathlib import Path
from fractions import Fraction
import json,re,collections

HERE=Path(__file__).resolve().parent
SPLITS=json.loads((HERE/"split_assignments.json").read_text())

def parse_utterance(path, speaker):
    text=path.read_text(errors='strict')
    items={}
    for line in text.split('Stream_Items\n',1)[1].split('End_of_Stream_Items',1)[0].splitlines():
        fields=line.split(';')
        first=fields[0].strip().split(maxsplit=1)
        if not first:continue
        i=int(first[0]);features={}
        fields[0]=first[1]
        for field in fields:
            pair=field.strip().split(maxsplit=1)
            if len(pair)==2:features[pair[0]]=pair[1].strip('"')
        items[i]=features
    relations={}
    for name,body in re.findall(r'Relation (\S+) ;[^\n]*\n(.*?)End_of_Relation',text,re.S):
        relations[name]={int(v[0]):list(map(int,v[1:])) for line in body.splitlines() if (v:=line.split())}
    def chain(rel,start):
        seen=set()
        while start:
            assert start not in seen
            seen.add(start);yield start
            start=rel[start][3]
    segrel=relations['Segment'];segments={};previous=Fraction(0)
    head=next(i for i,v in segrel.items() if v[4]==0)
    seq=[]
    for index,i in enumerate(chain(segrel,head)):
        item_id=segrel[i][0];item=items[item_id]
        end=Fraction(item['end'])*1000
        segments[item_id]={'start_ms':str(previous),'end_ms':str(end),'duration_ms':str(end-previous),'index':index,'phone':item['name']}
        seq.append((end,item['name']));previous=end
    lab=path.parents[2]/'lab'/(path.stem+'.lab')
    labels=[];lab_previous=Fraction(0)
    for line in lab.read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        t,_,phone=line.split();end=Fraction(t)*1000
        labels.append({'start_ms':str(lab_previous),'end_ms':str(end),'duration_ms':str(end-lab_previous),'index':len(labels),'phone':phone});lab_previous=end
    speech_labels=[x for x in labels if x['phone'] not in ('pau','sil')]
    speech_ids=[k for k,v in segments.items() if v['phone'] not in ('pau','sil')]
    if len(speech_labels)!=len(speech_ids):
        return [],collections.Counter({'utterance_phone_count_mismatch':1})
    if not all(segments[k]['phone']==label['phone'] and abs(Fraction(segments[k]['end_ms'])-Fraction(label['end_ms']))<=Fraction(1,100) for k,label in zip(speech_ids,speech_labels)):
        return [],collections.Counter({'utterance_phone_label_or_time_mismatch':1})
    # Use actual lab intervals; extra leading/trailing pauses need not be in
    # SylStructure. A pause inside a word leaves a noncontiguous lab index.
    for k,label in zip(speech_ids,speech_labels):segments[k]=label
    rel=relations['SylStructure'];words=[];exclusions=collections.Counter()
    roots=[i for i,v in rel.items() if v[1]==0 and v[4]==0 and items[v[0]].get('name')!='syl' and v[0] not in segments]
    assert len(roots)==1,(path,roots)
    for word_index,i in enumerate(chain(rel,roots[0])):
        v=rel[i];word=items[v[0]]['name']
        syls=[]
        for si in chain(rel,v[2]):
            sv=rel[si];feat=items[sv[0]]
            if feat.get('name')!='syl':raise ValueError((path,'non-syllable child'))
            phones=[segments[rel[pi][0]] for pi in chain(rel,sv[2])]
            syls.append({'stress':int(feat.get('stress',0)),'phones':[p['phone'] for p in phones],
                         'duration_ms':str(sum((Fraction(p['duration_ms']) for p in phones),Fraction(0))),
                         'segment_indices':[p['index'] for p in phones]})
        if len(syls)!=2:exclusions['not_two_syllables']+=1;continue
        if any(not s['phones'] for s in syls):exclusions['empty_syllable']+=1;continue
        inds=[j for s in syls for j in s['segment_indices']]
        if inds!=list(range(inds[0],inds[-1]+1)):exclusions['noncontiguous_phones']+=1;continue
        if any(p in ('pau','sil') for s in syls for p in s['phones']):exclusions['pause_inside_word']+=1;continue
        if any(Fraction(p['duration_ms'])<=0 for p in segments.values() if p['index'] in inds):exclusions['nonpositive_phone_duration']+=1;continue
        split=SPLITS[path.stem]
        words.append({'id':f'{speaker}/{path.stem}/{word_index}','speaker':speaker,'utterance':path.stem,'word_index':word_index,
                      'word':word,'split':split,'syllables':syls,'total_ms':str(sum((Fraction(s['duration_ms']) for s in syls),Fraction(0)))})
    return words,exclusions

def main():
    rows=[];audit={}
    for speaker in ['slt','bdl']:
        folder=HERE/'data'/f'cmu_us_{speaker}_arctic';count=collections.Counter();files=sorted((folder/'festival/utts').glob('*.utt'))
        for path in files:
            got,skips=parse_utterance(path,speaker);rows.extend(got);count.update(skips)
        sr=[x for x in rows if x['speaker']==speaker]
        audit[speaker]={'utterance_files':len(files),'retained_words':len(sr),'by_split':dict(collections.Counter(x['split'] for x in sr)),
                        'utterances_by_split':{s:len(set(x['utterance'] for x in sr if x['split']==s)) for s in ['train','validation','test']},'exclusions':dict(count)}
    (HERE/'prepared_words.json').write_text(json.dumps(rows,indent=2)+'\n')
    (HERE/'data_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2));print('First retained word',json.dumps(rows[0],indent=2))

if __name__=='__main__':main()
