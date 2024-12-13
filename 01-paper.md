---
title: Enriching paleogeographic reconstruction through the integration of Bayesian age models into site-based apparent polar wander paths
subtitle: Enriching APWPs
short_title: Bayesian age models for APWPs
authors:
  - name: Nicholas L. Swanson-Hysell
    affiliations:
      - Institute for Rock Magnetism, Department of Earth and Environmental Sciences, University of Minnesota
    orcid: 0000-0003-3215-4648
    email: nicks-h@umn.edu
  - name: Yiming Zhang
    affiliations:
      - Institute for Rock Magnetism, Department of Earth and Environmental Sciences, University of Minnesota
    orcid: 0000-0003-3215-4648
    email: nicks-h@umn.edu
license: CC-BY-4.0
keywords: paleogeography, polar wander, plate tectonics
---

**Abstract**

Apparent polar wander paths (APWPs) synthesized from paleomagnetic data are foundational to our understanding of the long-term motions of lithospheric plates. Recent progress has highlighted the advantages of developing such paths from paleomagnetic sites that provide individual snapshots of the ancient geomagnetic field. This approach contrasts with traditional methods that have relied on study level mean poles and through which it is challenging to incorporate positional and temporal uncertainty. A benefit of the site-based approach is the flexibility in incorporating temporal information. Rather than relying on a single age assigned to a formation mean paleomagnetic pole, varied ages can be assigned to individual sites. This flexibility enables the incorporation of temporal uncertainty. At the same time, it allows the incorporation of temporal information that would otherwise be discarded such as the constraint imposed by the principle of stratigraphic superposition --- a site stratigraphically above another is younger. In this contribution, we show how age models developed through Markov-Chain Monte Carlo methods can be incorporated into 

## Background

Reconstructing the movement of Earth's lithospheric plates through time is central to our understanding of the long-term dynamics of Earth's surface and interior. Paleomagnetic data are central to our understanding of plate motions throughout Earth history. A central aspect of synthesizing paleomagnetic data for paleogeographic reconstruction is the construction of apparent polar wander paths (APWPs). Apparent polar wander refers to the appearance in paleomagnetic datasets that the position of the pole has moved relative to the continent. In actuality, the pole stayed in the same position while the plate moved relative to the pole. Reconstructing the path of the apparent pole motion through an APWP provides the basis for reconstructing plate motion relative to the pole.

Historically, paleomagnetists have sought to develop mean paleomagnetic poles associated with geologic formations. These poles comprise a collection of individual paleomagnetic sites with the goal of averaging secular variation of Earth's geomagnetic field. Paleomagnetic directions from individual sites can be transformed into ``virtual geomagnetic poles'' (VGPs) that represent the position of geomagnetic north at the time that site magnetization was acquired. Using the time-averaged geocentric axial dipole hypothesis, the mean of these VGPs is taken to correspond to the spin axis (either the geographic north or geographic south pole). Paleomagnetic poles are reported with a 95% confidence bound on the mean position calculated using Fisher statistics. For this mean paleomagnetic pole to be an accurate representation of geographic north, it is necessary for sufficient number of site data to have been incorporated to average out secular variation of the geomagnetic field about the spin axis. Various criteria have been introduced to evaluate whether this goal has been achieved ([](doi:10.1111/j.1365-246X.2011.05050.x)).

Recent work has highlighted how the step of grouping sites into means prior to APWP synthesis has multiple deficiences particularly in giving equal weight to mean poles that can themselves be informed by vastly different numbers of sites ([](10.1029/2023GL103436)). Subsequent work has revealed how constructing apparent polar wander path from site level data

The []10.17605/osf.io/TQX3F

Scientific communication today is designed around print documents and pay-walled access to content. Over the last decade, the open-science movement has accelerated the use of pre-print services and data archives that are vastly improving the accessibility of scientific content. However, these systems are not designed for communicating modern scientific outputs, which encompasses **so much more** than a paper-centric model of the scholarly literature.

> We believe how we share and communicate scientific knowledge should evolve past the status quo of print-based publishing and all the limitations of paper.

The communication and collaboration tools that we are building in the Project Jupyter are built to follow the FORCE11 recommendations (Bourne _et al._, 2012). Specifically:

1. rethink the unit and form of scholarly publication;
2. develop tools and technologies to better support the scholarly lifecycle; and
3. add data, software, and workflows as first-class research objects.

By bringing professional, high-quality tools for science communication into the research lifecycle, we believe we can improve the collection and preservation of scholarly metadata (citations, cross-references, annotations, etc.) as well as open up new ways to communicate science with interactive figures & equations, computation, and reactivity.

The tools that are being built by the Project Jupyter are focused on introducing a new Markup language, MyST (Markedly Structured Text), that works seamlessly with the Jupyter community to enhance and promote a new path to document creation and publishing for next-generation scientific textbooks, blogs, and lectures. Our team is currently supported by the [Sloan Foundation](https://sloan.org), ([Grant #9231](https://sloan.org/grant-detail/9231)).

MyST enables rich content generation and is a powerful format for scientific and technical communication. JupyterBook uses MyST and has broad adoption in publishing tutorials and educational content focused around Jupyter Notebooks.

> The components behind Jupyter Book are downloaded 30,000 times a day, with 750K downloads last month.

The current toolchain used by [JupyterBook] is based on [Sphinx], which is an open-source documentation system used in many software projects, especially in the Python ecosystem. `mystjs` is a similar tool to [Sphinx], however, designed specifically for scientific communication. In addition to building websites, `mystjs` can also help you create scientific PDFs, Microsoft Word documents, and JATS XML (used in scientific publishing).

`mystjs` uses existing, modern web-frameworks in place of the [Sphinx] build system. These tools come out-of-the-box with prefetching for faster navigation, smaller network payloads through modern web-bundlers, image optimization, partial-page refresh through single-page application. Many of these features, performance and accessibility improvements are difficult, if not impossible, to create inside of the [Sphinx] build system.

In 2022, the Executable Books team started work to document the specification behind the markup language, called [myst-spec](https://github.com/jupyter-book/myst-spec), this work has enabled other tools and implementations in the scientific ecosystem to build on MyST (e.g. [scientific authoring tools](https://curvenote.com/for/writing), and [documentation systems](https://blog.readthedocs.com/jupyter-book-read-the-docs/)).

The `mystjs` ecosystem was developed as a collaboration between [Curvenote], [2i2c] and the [ExecutableBooks] team. The initial version of `mystjs` was originally release by [Curvenote] as the [Curvenote CLI](https://curvenote.com/docs/cli) under the MIT license, and transferred to the [ExecutableBooks] team in October 2022. The goal of the project is to enable the same rich content and authoring experiences that [Sphinx] allows for software documentation, with a focus on web-first technologies (Javascript), interactivity, accessibility, scientific references (e.g. DOIs and other persistent IDs), professional PDF outputs, and JATS XML documents for scientific archiving.

## MyST Project

In this paper we introduce `mystjs`, which allows the popular MyST Markdown syntax to be run directly in a web browser, opening up new workflows for components to be used in web-based editors, [directly in Jupyter](https://github.com/jupyter-book/jupyterlab-myst) and in JupyterLite. The libraries work with current MyST Markdown documents/projects and can export to [LaTeX/PDF](https://myst.tools/docs/mystjs/creating-pdf-documents), [Microsoft Word](https://myst.tools/docs/mystjs/creating-word-documents) and [JATS](https://myst.tools/docs/mystjs/creating-jats-xml) as well as multiple website templates using a [modern](https://myst.tools/docs/mystjs/accessibility-and-performance) React-based renderer. There are currently over 400 scientific journals that are supported through [templates](https://github.com/myst-templates), with [new LaTeX templates](https://myst.tools/docs/jtex/create-a-latex-template) that can be added easily using a Jinja-based templating package, called [jtex](https://myst.tools/docs/jtex).

In our paper we will give an overview of the MyST ecosystem, how to use MyST tools in conjunction with existing Jupyter Notebooks, markdown documents, and JupyterBooks to create professional PDFs and interactive websites, books, blogs and scientific articles. We give special attention to the additions around structured data, standards in publishing (e.g. efforts in representing Notebooks as JATS XML), rich [frontmatter](https://myst.tools/docs/mystjs/frontmatter) and bringing [cross-references](https://myst.tools/docs/mystjs/cross-references) and [persistent IDs](https://myst.tools/docs/mystjs/external-references) to life with interactive hover-tooltips ([ORCID, RoR](https://myst.tools/docs/mystjs/frontmatter), [RRIDs](https://myst.tools/docs/mystjs/external-references#research-resource-identifiers), [DOIs](https://myst.tools/docs/mystjs/citations), [intersphinx](https://myst.tools/docs/mystjs/external-references#intersphinx), [wikipedia](https://myst.tools/docs/mystjs/external-references#wikipedia-links), [JATS](https://myst.tools/docs/mystjs/typography), [GitHub code](https://myst.tools/docs/mystjs/external-references#github-links), and more!). This rich metadata and structured content can be used directly to improve science communication both through self-publishing books, blogs, and lab websites — as well as journals that incorporate Jupyter Notebooks.

## Features of MyST

MyST is focused on scientific writing, and ensuring that citations are first class both for writing and for reading (see Figure 1).

:::{figure} ./images/citations.png
:label: citations
Citations are rendered with a popup directly inline.
:::

MyST aims to show as much information in context as possible, for example, Figure 2 shows a reading experience for a referenced equation: you can immediately **click on the reference**, see the equation, all without loosing any context -- ultimately saving you time. Head _et al._ (2021) found that these ideas both improved the overall reading experience of articles as well as allowed researchers to answer questions about an article **26% faster** when compared to a traditional PDF!

![](./images/equations.gif)
**Figure 2**: In context cross-references improve the reading experience.

One of the important underlying goals of practicing reproducibility, sharing more of the methods and data behind a scientific work so that other researchers can both verify as well as build upon your findings. One of the exciting ways to pull for reproducibility is to make documents directly linked to data and computation! In Figure 3, we are showing outputs from a Jupyter Notebook directly part of the published scientific narrative.

![](./images/interactive.gif)
**Figure 3**: Embedding data, interactivity and computation into a MyST article.

To drive all of these features, the contents of a MyST document needs to be well defined. This is critical for powering interactive hovers, linked citations, and compatibility with scientific publishing standards like the Journal Article Metadata Tag Suite (JATS). We have an emerging specification for MyST, [`myst-spec`](https://spec.myst.tools), that aims to capture this information and transform it between many different formats, like PDF, Word, JSON, and JATS XML (Figure 4). This specification is arrived at through a community-centric MyST Enhancement Proposal ([MEP](https://compass.executablebooks.org/en/latest/meps.html)) process.

![](./images/structured-data.gif)
**Figure 4**: The data behind MyST is **structured**, which means we can transform it into many different document types and use it to power all sorts of exciting features!

One of the common forms of scientific communication today is through PDF documents. MyST has excellent support for creating PDF documents, using a data-driven templating library called `jtex`. The document in Figure 5 was created using MyST!

![](./images/pdf-two-column.png)
**Figure 5**: A PDF rendering through MyST.

## Conclusion

There are many opportunities to improve open-science communication, to make it more interactive, accessible, more reproducible, and both produce and use structured data throughout the research-writing process. The `mystjs` ecosystem of tools is designed with structured data at its core. We would love if you gave it a try -- learn to get started at <https://myst.tools>.
