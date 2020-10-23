---
comments: true
title: Power analysis for Fisher and Z methods
date: '2020-10-22 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Power
  - Fisher
  - Z
categories:
  - WGS data analysis
---

In this post I performed the power analysis for Fisher and weighted Z-test (also called ‘Stouffer’s method, see [Whitlock (2005)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x)) approaches. The major difference between these two methods is that during the p-value combination step, the Z method favours symmetric rejection and is less sensitive to a single low p-value, requiring more consistently low p-values to yield a low combined p-value.

Key features in the simulation,

- Delta_p Levene's model




