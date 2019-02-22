reserved_atoms = {'Br','Cl','C','N','O','S','P','F','I','B'}
unresered_atom_mark = 'X'
    
def smarts2words(sm): 
    qmol = Chem.MolFromSmarts(sm)
    l,mtos = atomnfp(qmol)
    qmol = Chem.MolFromSmiles(mtos)
    s = Chem.MolToSmarts(qmol)
    k_start =[]
    k_end = []
    for n,c in enumerate(s):
        if c == '[':
            k_start.append(n)
        elif c == ']':
            k_end.append(n)
    r=[]
    count = 0
    r.append(s[0:k_start[0]])
    r.append(l[count][1])
    count += 1
    for n in range(len(k_start)-1):
        r.append(s[k_end[n]+1:k_start[n+1]])
        r.append(l[count][1])
        count += 1
    r.append(s[k_end[-1]+1:])
    return ' '.join(r).strip()

def atomnfp(qmol):
    Chem.Cleanup(qmol)
    Chem.SanitizeMol(qmol)
    qmol.UpdatePropertyCache()
    amp = [atom.GetAtomMapNum() for atom in qmol.GetAtoms()]
    for atom in qmol.GetAtoms():atom.SetAtomMapNum(0)
    mtos = Chem.MolToSmiles(qmol)
    sao = qmol.GetProp('_smilesAtomOutputOrder')[1:-2].split(',')
    r =[]
    for n,i in enumerate(sao):
        atom = qmol.GetAtomWithIdx(int((i.strip())))        
        nfp=[]
        nb = []
        if amp[n] != 0:amp[n] = 1
        nfp.append(str(amp[n]))
        sb = atom.GetSymbol()
        if sb in reserved_atoms:
            nfp.append(sb)
        else:
            nfp.append(unresered_atom_mark)           
        nfp.append('$')
        nfp.append(str(atom.GetTotalValence()))
        for n in atom.GetNeighbors():
            sb = atom.GetSymbol()
            if sb in reserved_atoms:
                nb.append(sb)
            else:
                nb.append(unresered_atom_mark) 
        nfp.extend(sorted(nb))
        nfp.append(str(atom.GetTotalNumHs()))
        r.append([atom.GetSymbol(),''.join(nfp)])
    return r,mtos
