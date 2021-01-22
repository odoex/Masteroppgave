# Masteroppgave

## Advection equation
Den stabile mesh refinement metoden fungerer nå for Advection equation. (Tror jeg)
Den ligger under Advection_equation\stableMethod. Tilhørende filer er: 
- meshRefinementAdv2.m
- rhsFineGridAdv.m
- projectFineBoundary.m

Tilhørende filer som også brukes i vanlig mesh refinement advection er:
- createSolutionVector.m
- finiteVolumeAdv.m
- RK_4Adv.m
- exactSolAdv.m
- f_adv.m
- g_adv.m
- rhsAdv.m
- boundary.m

De to førstenevnte, finiteVolumeAdv.m og RK_4Adv.m,
må modereres for å bruke den stabile metoden. For finiteVolumeAdv.m 
byttes den vanlige projectFlux() ut med projectFineBoundary() (rart navn, 
skal byttes). For RK_4.m kommenteres den vanlige rk-metoden ut og den 
nedenforstående metoden brukes, hvor det sjekkes om gridet er et undergrid
hvor det i så fall kalles på funksjonen rhsFineGridAdv() som håndterer 
overgangen mellom fint og grovt grid på en annen måte enn som en rand. To 
av de tre testene i mappen tests fungerer ikke selv om det er riktig. Må
finne ut av hvorfor.

## Euler equations
### Normal mesh refinement
Euler ligningene sin mesh refinement metode fungerer så vidt jeg vet, men 
gir ikke ønskede resultater ettersom den har  mye lavere konvergens. Den 
befinner seg i mappen Euler_equations. Tilhørende filer er: 
- meshRefinement.m
- exactSolEuler.m
- finiteVolume.m
- rhs.m
- boundary.m
- f.m
- g.m
For konstant løsning (disse fungerer i stede for exactSolEuler.m:
- boundaryValuesEulerConstant.m
- initialConditionsEulerConstant.m
Fra mappen meshRefinement som inneholder filer som kan brukes av alle 
metodene er også: 
- Node.m
- projectFlux.m
- RK_4.m

Metoden fungerer på følgende måte: Som vanlig skapes hovedgrid og undergrid
i meshRefinement filen. exactSolEuler.m gir eksakt løsning til gridene ved 
at to vektorer x og y sendes inn i funksjonen, som spanner gridet som skal 
fylles med eksakte verdier. Så starter metoden i finiteVolume.m. Her 
opdaterer metoden nettet til nåverende hoved-node, som til å starte med er 
det grove hovednettet, og lagrer det i variabelen u. (Dette er strengt tatt
ikke nødvendig før etterpå, men ble gjort slik fordi jeg leste en oppskrift 
i denne rekkefølgen.) Videre sjekkes det om G har et undergrid, og i så 
fall kjøres finite volumes metoden rekursivt om igjen med undergridet 
G.child som hovedgrid. Slik fortsetter det til et av gridene ikke har noen 
undergrid. Da settes gridet i noden G.u til å være lik u, slik at det er 
opdatert et tidssteg. Hvis dette tidssteget er mindre en foreldregridet 
sitt tidssteg, loopes det gjennom tidsstegene og opdaterer ved hvert 
tidssteg helt til tiden er ved foreldregrides neste tidssteg, der 
foreldregridet sist ble opdatert (selv om ikke selve noden ble opdatert 
enda, det skjer når alle barna er opdatert slik at barna skal opdateres fra
forrige tidssteg og ikke neste). Deretter hopper metoden tilbake til 
forelderen og opdaterer dette gridet i selve noden. Etter hver gang et av 
foreldregridene er opdatert kalles funksjonen projectFlux.m hvor alle
verdiene fra det fine gridet som overlapper punkter på det grovere 
foreldregridet overføres til foreldregridet. 

Når nettet opdateres sendes det inn i RK_4 som er Runge Kutta metoden. Den 
bruker igjen rhs.m som er finite volumes metoden. Denne bruker igjen 
boundary.m metoden, som henter ekte randdata fra exactSolEuler.m, og 
interpolerer for å få de manglende verdiene på randen når det er 
forfiningen som opdateres. f.m og g.m fungerer som flux funksjoner slik at
under finite volumes utregningene transformeres alle verdiene i gridet u 
til de tilhørende flux verdiene i x-retning, f, og y-retning, g. 

### Stable mesh refinement:
For stabil mesh refinement av Euler equations metoden er følgende metoder 
opprettet: 
- stableMeshRefinement.m
- initiateSubgrid.m
- finiteVolumeStableMethod.m
- RK_4StableMethod.m
- rhsFineGrid.m
- diffusion.m
- projectFluxStableMethod.m
In addition to the following files from the Euler equations folder:
- exactSolEuler.m
- rhs.m
- f.m
- g.m
- boundary.m
and from meshRefinement folder:
- Node.m

Stable mesh refinement fungerer som vanlig mesh refinement med følgende 
unntak: 

stableMeshRefinement.m er laget av meshRefinement.m med moderasjonene til 
meshRefinementAdv2.m, som er den stabile mesh refinement metoden for 
advection equation, slik at det er tilrettelagt for et ekstra lag punkter i 
det fine gridet som trengs for denne framgangsmåten. initiateSubgrid.m 
fungerer nå og gir initialdata til alle undergrid, og legger på de ekstra 
punktene som kun subgrids skal ha. Den benytter seg også av exactSolEuler.m
for initialverdiene. Tilhørende test initiateSubgridTest.m i mappen tests
fungerer også. 

I RK_4StableMethod.m deles det inn i to tilfeller: For hovedgridet tas det 
ikke hensyn til andre nett og denne sendes bare inn i den vanlige rhs.m 
metoden. Alle undergrid sendes inn i rhsFineGrid.m hvor det blir nøye fyllt
inn i punktene utenfor gridet, på boundaryen og inni gridet slik at 
volumene er laget på en måte som gjør metoden stabil. I tillegg kalles det 
på funksjonen diffusion.m som tar i mot et nett u og regner ut diffusjonen 
som skal legges til flux-utregningene i dette punktet. Dette gjøres for 
alle punkter inni metoden, og det returneres en matrise U_diff med disse 
verdiene som legges til fluxen. Diffusjonsfunksjonen regner ut to matriser 
U_pluss og U_minus som tar for seg den inngående diffusjonen i de to 
fluxene f_i+1/2 og f_i-1/2 som til slutt legges sammen og sendes tilbake i 
U_diff. 
