# Standortanalyse Rebbau


This Repo holds the code to the project [Standortanalyse Rebbau in den
Kantonen Aargau, Bern und
Zug](https://www.zhaw.ch/en/research/project/74192) (link not active
yet).

It holds several scripts:

1.  `get_swissalti3d.R`: To get swissalti3D data, we obtained a csv with
    URLs to the download path (`data/swissalti3d_2m_all.csv`). This
    Scripts downloads all (?) files to `data/swissAlti3D/`
2.  `historic_HI.R`: Calculates the historic (measured) *Hugglin Index*
    for the three cantons based on the `*.nc` files in
    `data/Klimadaten_Feb24`.
3.  `future_HI.R`: Calculates the expected *Hugglin Index* in the future
    based on modeled data
4.  `HI-pro-rebberg.R`: Extract the *Hugglin Index* per vineyard (not in
    production)

See also this previous project [Standortanalyse von
pilzwiderstandsfähigen Rebsorten für den Kanton
Luzern](https://digitalcollection.zhaw.ch/items/8a2102f6-582f-482b-8301-64cd60fbf45f).
