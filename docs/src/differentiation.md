```@meta
CurrentModule = NaturalNeighbours
```

# Differentiation Example 

The purpose of this example is to explore derivative generation. For this, it is important to note that we are thinking of _generating_ derivatives rather than _estimating_ them: Following [Alfeld (1989)](https://doi.org/10.1016/B978-0-12-460515-2.50005-6), derivative generation only seeks to find derivatives that best fit our assumptions of the data, i.e. that give a most satisfactory interpolant, rather than trying to find exact derivative values. The complete quote for this by [Alfeld (1989)](https://doi.org/10.1016/B978-0-12-460515-2.50005-6) is below:

> It seems inevitable that in order to obtain an interpolant that is both local and smooth one has to supply derivative data. Typically, such data are not part of the interpolation problem and have to be made up from existing func tional data. This process is usually referred as derivative estimation, but this is probably a misnomer. The objective is not to estimate existing but unknown values of derivatives. Instead, it is to generate values that will yield a satisfactory interpolant. Even if an underlying primitive function did exist it might be preferable to use derivative values that differ from the exact ones. (For example, a maximum error might be decreased by using the "wrong" derivative values.) Therefore, I prefer the term derivative generation rather than derivative estimation.