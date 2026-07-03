# AFM

Progetto di analisi dati AFM (Atomic Force Microscopy) — Università di Milano Bicocca, Monza.
Obiettivo scientifico: misurare le proprietà meccaniche (rigidezza/modulo elastico) di cellule
tumorali a partire da curve forza-distanza raccolte con un cantilever AFM, confrontando linee
cellulari e regioni di misura diverse.

## Struttura del repository

```
Code_General/
    01_JPK_reader_and_preliminary_QC.ipynb   # notebook attuale, "definitivo" per lo stage 1
    Preliminary/                              # versioni esplorative/precedenti (v0..v4), superate
    Curve di forza belle/                     # dati grezzi .jpk-force (input), organizzati a gruppi
    afm_output/                               # output generato dal notebook 01 (git-ignored di fatto)
Docs/
    AFM_ML.pdf, AFM_presentation.pdf          # slide/relazioni del progetto
    afm_ml_stat_project_presentation.html     # presentazione standalone: proof of concept ML su dati sintetici
    pipeline_qc_afm_documento_guida.pdf       # guida di riferimento per la pipeline di QC
```

Non ci sono script di build/test: il progetto è una pipeline di notebook Jupyter eseguiti in ordine,
tipico lavoro di analisi dati scientifica in Python (numpy/pandas/matplotlib).

## Dati grezzi

`Code_General/Curve di forza belle/<gruppo>/*.jpk-force` — file nativi JPK (in realtà zip contenenti
`header.properties` + segmenti). Quattro gruppi sperimentali attualmente presenti:

- `Petri20k_MCF7R_t48h_cito_20260225` (55 file)
- `Petri20k_MCF7R_t48h_nucleo_20260225` (60 file)
- `T47DR_t48h_cito_20260318` (59 file)
- `T47DR_t48h_nucleo_20260318` (83 file)

Il nome del gruppo codifica: linea cellulare (`MCF7R` / `T47DR`, entrambe linee di carcinoma
mammario), regione di misura (`cito` = citoplasma, `nucleo` = nucleo), tempo (`t48h`) e data
dell'esperimento (`YYYYMMDD`). Il parsing di questi metadati dal nome cartella avviene in
`parse_group_name()` nel notebook 01.

Ogni file `.jpk-force` contiene almeno 2 segmenti: `segment 0` = extend/approach (candidato per il
fit meccanico), `segment 1` = retract (utile per adesione/isteresi).

## Pipeline di analisi (stage previsti)

1. **01_JPK_reader_and_preliminary_QC.ipynb** (fatto, in `Code_General/`)
   Legge tutti i `.jpk-force`, esporta ogni segmento in CSV, salva metadata/proprietà JPK,
   costruisce il summary delle curve `segment0/extend` e assegna un **QC preliminare
   strumentale** (z-score robusti globali e per gruppo su feature come rumore di baseline,
   escursione di forza, salti di segnale). Importante: **nessuna curva viene scartata in questo
   notebook** — il QC qui è solo diagnostico/informativo, serve a fondare le decisioni successive.
   Produce `afm_output/manifests/segment0_manifest.csv`, che è l'input del notebook successivo.

2. **02_preprocessing_contact_point.ipynb** (non ancora scritto in `Code_General/`; esiste solo
   una serie di prototipi in `Preliminary/02_contact_point_and_hertz_like_fit_v1..v4.ipynb`)
   Dovrà occuparsi di: refinement della baseline, stima del contact point, costruzione
   dell'indentazione, QC sul contact point, fit meccanico Hertz-like `F(δ) = A·δ^(3/2)` (A =
   parametro elastico efficace, non ancora modulo di Young assoluto — servirebbe raggio della
   punta e Poisson ratio per la conversione), curve forza-indentazione preprocessate e plot
   diagnostici curva-per-curva.

3. **03 — statistica dei fit** (solo prototipo in `Preliminary/03_fit_results_statistics.ipynb`)
   Analizza la distribuzione di A (tipicamente lavorando su `log A`, più stabile), confronta i
   gruppi sperimentali, stima intervalli bootstrap sulle mediane di gruppo, a partire solo dai fit
   di qualità "good" (le curve "suspicious" sono escludibili/includibili a scelta).

La cartella `Preliminary/` contiene le versioni esplorative superate dei notebook 02/03: utile
come riferimento sulla logica già sviluppata, ma non è il codice "vivo" — quando si riscrive 02 o
03 in `Code_General/`, ripartire dalla logica lì ma ripulirla/consolidarla come è stato fatto per 01.

## Convenzioni del codice

- Notebook e commenti sono in **italiano**; mantenere questa lingua per coerenza.
- Ogni notebook definitivo inizia con una cella markdown che dichiara esplicitamente cosa fa e
  cosa **non** fa (v. cella 0 di 01), e termina con una cella che descrive l'output prodotto e il
  notebook successivo atteso — mantenere questo pattern nei nuovi notebook.
- Il notebook 01 scrive un `afm_output/README_outputs.txt` che documenta ogni sottocartella di
  output generata (`csv_curves/`, `metadata/`, `summaries/`, `manifests/`, `plots_overview/`,
  `plots_qc/`) — aggiornarlo se cambia lo schema di output.
- QC preliminare (notebook 01) è filosoficamente diverso dal QC del notebook 02/03: il primo è
  solo diagnostico e non esclude nulla; l'esclusione delle curve va fatta più a valle, dopo il fit,
  in modo esplicito e tracciabile.
- `afm_output/` è interamente rigenerabile dal notebook 01 a partire da `Curve di forza belle/`:
  trattarla come cache, non come sorgente di verità.

## Componente ML/statistica sui dati sintetici

`Docs/afm_ml_stat_project_presentation.html` documenta un mini-progetto propedeutico, separato dai
notebook di `Code_General/`, su **1000 curve AFM sintetiche** (950 buone + 50 con difetti
strumentali iniettati ad arte: `no_contact`, `saturation`, `huge_noise`, `jump`, `wrong_slope`).
Serve come dimostrazione/validazione della pipeline generale — curve grezze → correzione baseline →
stima contatto → feature fisiche → anomaly detection → fit meccanico → statistica dei parametri —
prima di applicarla ai dati reali JPK. Utile come riferimento concettuale quando si progetta il QC
o l'anomaly detection reale, ma il codice lì non è collegato direttamente ai notebook di
`Code_General/`.
