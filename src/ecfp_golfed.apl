ECFP←{
 ⍺←2 ⋄ r←⊃⍺ ⋄ S←2⊃⍺,1 ⋄ w←⎕CSV ⍵ ⍬ 4 ⋄ n←⊃⌽⍴w
 P←S×n↑((≢w)>n)/w[≢w;] ⋄ M←n↑w ⋄ D←1 1⍉M ⋄ B←M×∘.≠⍨⍳n ⋄ A←0≠B
 K←⍸A∧(⍳n)∘.<⍳n ⋄ I←(⍳n)∘.∊K ⋄ m←≢K
 F←{
  (t d q c aI aC)←⍵
  N←{⍸0≠B[⍵;]}¨⍳n
  s←S∧(P≠0)∧(~q)∧{(≢⍵)=≢∪d[⍵]}¨N
  g←{e←⍵⊃N ⋄ p←⍉↑(B[⍵;e])(d[e]) ⋄ {2147483647|⍺+31×⍵}/⌽1,(t,(d[⍵]),,p[⍋p;]),s[⍵]/P[⍵]}¨⍳n
  y←c∨I∨A∨.∧c
  (t+1)g(q∨s)y(aI,g)(aC⍪y)
 }
 z←(F⍣r)1 D(n⍴0)((n,m)⍴0)D((n,m)⍴0)
 aI←5⊃z ⋄ aC←6⊃z ⋄ e←0=+/aC ⋄ b←∪e/aI
 T←(~e)/n/0,⍳r ⋄ aI←(~e)/aI ⋄ aC←(~e)⌿aC
 o←⍋⍉↑T aI ⋄ aI←aI[o] ⋄ aC←aC[o;]
 {⍵[⍋⍵]}∪b,aI/⍨(⍳≢aC)=(↓aC)⍳↓aC
}