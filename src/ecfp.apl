Hash←{{2147483647|⍺+31×⍵}/⌽1,⍵}
ReadMatrix ← {⎕CSV ⍵ ⍬ 4}

ECFP←{
    ⍺←2                                              ⍝ default radius, stereo on
    radius←⊃⍺
    useStereo←2⊃⍺,1
    raw←ReadMatrix ⍵
    n←⊢/⍴raw                                         ⍝ atoms = number of columns
    parity←useStereo×(1+(≢raw)>n)⊃(n⍴0)(raw[≢raw;])  ⍝ optional last row
    M←n↑raw                                          ⍝ the square matrix itself
    ids0←1 1⍉M                                       ⍝ diagonal → initial identifiers
    B←M×∘.≠⍨⍳n                                       ⍝ zero the diagonal → bond orders
    A0←0≠B                                           ⍝ boolean adjacency
    BP←⍸A0∧(⍳n)∘.<⍳n                                 ⍝ upper-triangle bonds as (i j)
    m←≢BP
    incident←(⍳n)∘.{⍺∊⍵}BP                           ⍝ n×m : does atom a touch bond k?

    Step←{                                           ⍝ one iteration, all atoms at once
        (t ids done cover accID accIT accCOV)←⍵
        distinct←{nb←⍸0≠B[⍵;] ⋄ (≢nb)=≢∪ids[nb]}¨⍳n
        should←useStereo∧(parity≠0)∧(~done)∧distinct
        newids←{nb←⍸0≠B[⍵;] ⋄ P←⍉↑(B[⍵;nb])(ids[nb]) ⋄ seq←t,(ids[⍵]),,P[⍋P;] ⋄ Hash seq,should[⍵]/parity[⍵]}¨⍳n
        newcover←cover∨incident∨A0∨.∧cover
        (t+1)newids(done∨should)newcover(accID,newids)(accIT,n⍴t)(accCOV⍪newcover)
    }

    init←1 ids0(n⍴0)((n,m)⍴0)ids0(n⍴0)((n,m)⍴0)
    final←(Step⍣radius)init
    accID←5⊃final ⋄ accIT←6⊃final ⋄ accCOV←7⊃final

    empt←0=+/accCOV                                  ⍝ empty cover ≡ radius-0 feature
    base←∪empt/accID                                 ⍝ radius-0: dedup by identifier
    keep←(~empt)/⍳≢accID
    ID2←accID[keep] ⋄ IT2←accIT[keep] ⋄ COV2←accCOV[keep;]
    ord←⍋⍉↑IT2 ID2                                   ⍝ order by (iteration, identifier)
    ID2←ID2[ord] ⋄ COV2←COV2[ord;]
    first←(↓COV2)⍳↓COV2                              ⍝ first feature sharing each cover
    win←first=⍳≢first                                ⍝ keep the min (iteration, id)
    {⍵[⍋⍵]}∪base,win/ID2                             ⍝ unique identifiers, sorted
}
