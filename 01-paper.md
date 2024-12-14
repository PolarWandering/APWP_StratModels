---
title: Enriching paleogeographic reconstruction through the integration of Bayesian age models into site-based apparent polar wander paths
subtitle: Enriching APWPs
short_title: Bayesian age models for APWPs
authors:
  - name: Yiming Zhang
    affiliations:
      - Institute for Rock Magnetism, Department of Earth and Environmental Sciences, University of Minnesota
    orcid: 0000-0002-1407-302X
    email: yiming-z@umn.edu
  - name: Nicholas L. Swanson-Hysell
    affiliations:
      - Institute for Rock Magnetism, Department of Earth and Environmental Sciences, University of Minnesota
    orcid: 0000-0003-3215-4648
    email: nicks-h@umn.edu
  - name: Facu Sapienza
    affiliations:
      - Department of Geophysics, Stanford University
    orcid: 0000-0003-4252-7161
    email: sapienza@stanford.edu
  - name: Mat Domeier
    affiliations:
      - Centre for Planetary Habitability, University of Oslo
    orcid: 0000-0002-7647-6852
    email: m.m.domeier@geo.uio.no
  - name: Leandro Gallo
    affiliations:
      - Centre for Planetary Habitability, University of Oslo
    orcid: 0000-0002-8124-7536
    email: l.c.gallo@geo.uio.no

license: CC-BY-4.0
keywords: paleogeography, polar wander, plate tectonics
---

**Abstract**

Apparent polar wander paths (APWPs) synthesized from paleomagnetic data are foundational to our understanding of the long-term motions of lithospheric plates. Recent progress has highlighted the advantages of developing such paths from paleomagnetic sites that provide individual snapshots of the ancient geomagnetic field. This approach contrasts with traditional methods that have relied on study level mean poles and through which it is challenging to incorporate positional and temporal uncertainty. A benefit of the site-based approach is the flexibility in incorporating temporal information. Rather than relying on a single age assigned to a mean paleomagnetic pole, varied ages can be assigned to individual sites. This flexibility enables the incorporation of temporal uncertainty. At the same time, it allows the incorporation of temporal information that would otherwise be discarded such as the constraints imposed by the principle of stratigraphic superposition — a site stratigraphically above another is younger. In this contribution, we show how stratigraphic age models developed through Markov-Chain Monte Carlo methods can be integrated into the development of the site-based APWP. Through resampling paths from the posterior distributions of such age models, the temporal uncertainty can be incorporated while honoring stratigraphic superposition. We illustrate this approach using volcanostratigraphic successions with U-Pb age constraints from the late Mesoproterozoic Midcontinent Rift of North America that we use to develop a stratigraphically constrained site-based APWP of the late Mesoproterozoic Keweenawan Track.

## Background

Reconstructing the movement of Earth's lithospheric plates through time is central to our understanding of the long-term dynamics of Earth's surface and interior. Paleomagnetic data are central to our understanding of plate motions throughout Earth history. A central aspect of synthesizing paleomagnetic data for paleogeographic reconstruction is the construction of apparent polar wander paths (APWPs). Apparent polar wander refers to the appearance in paleomagnetic datasets that the position of the pole has moved relative to the continent. In actuality, the pole stayed in the same position while the plate moved relative to the pole. Reconstructing the path of the apparent pole motion through an APWP provides the basis for reconstructing plate motion relative to the pole.

Historically, paleomagnetists have sought to develop mean paleomagnetic poles associated with geologic formations. These poles comprise a collection of individual paleomagnetic sites with the goal of averaging secular variation of Earth's geomagnetic field. Paleomagnetic directions from individual sites can be transformed into ``virtual geomagnetic poles'' (VGPs) that represent the position of geomagnetic north at the time that site magnetization was acquired. Using the time-averaged geocentric axial dipole hypothesis, the mean of these VGPs is taken to correspond to the spin axis (either the geographic north or geographic south pole). Paleomagnetic poles are reported with a 95% confidence bound on the mean position calculated using Fisher statistics. For this mean paleomagnetic pole to be an accurate representation of geographic north, it is necessary for sufficient number of site data to have been incorporated to average out secular variation of the geomagnetic field about the spin axis. Various criteria have been introduced to evaluate whether this goal has been achieved ([](doi:10.1111/j.1365-246X.2011.05050.x)).

Recent work has highlighted how the step of grouping sites into means prior to APWP synthesis has multiple deficiences particularly in giving equal weight to mean poles that can themselves be informed by vastly different numbers of sites ([](10.1029/2022JB023953); [](10.1029/2023GL103436)).

The []10.17605/osf.io/TQX3F


:::{figure} ./images/citations.png
:label: citations
Citations are rendered with a popup directly inline.
:::
