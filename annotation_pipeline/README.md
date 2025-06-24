# IMPC Annotation Pipeline
The IMPC annotation pipeline assigns Mammalian Phenotype (MP) terms to significant genetic effects based on a p-value threshold of 0.0001. The goal is to associate phenotypic observations with the corresponding genetic modifications.
At the IMPC, genetic effects are identified using three statistical analysis platforms:
1. Linear Mixed Model framework or MM
2. Fisher's Exact Test framework or FE
3. Reference Range Plus Test framework or RR

## Continuous data
Continuous data are typically analysed using a linear mixed model framework. These continuous measurements are particularly informative because the direction of change can be determined through the effect size.

However, due to the complexity of the data, not all continuous variables can be analysed using this framework. In such cases, the IMPC often employs the Reference Range Plus (RR) method. Control data are first discretised into three categories: low, normal, and high. Mutant data points are then classified into one of these reference categories. Finally, a Fisher's Exact Test is applied to determine whether there is a statistically significant deviation from the normal category.

## Categorical data 
Categorical data in the IMPC encompasses a range of qualitative measurements and is analysed using Fisher’s Exact Test, as implemented in the R package OpenStats.

# How IMPC Annotation Pipeline Works
The `annotationChooser` function processes statistical analysis results called statpackets. It determines calls based on significance levels. These calls are then mapped to Mammalian Phenotype (MP) ontology terms using a provided `mp_chooser_file`. Finally, it updates the input statpacket's JSON component with the identified MP terms. If no relevant annotation is found or the statistical result is not significant, it returns the original statpacket with no MP terms added.

The annotation pipeline requires a reference table that summarises the available MP terms for a given IMPC parameter. This reference can be retrieved from [IMPReSS](https://www.mousephenotype.org/impress/index).
The ETL pipeline handles this by generating the `mp_chooser.json` file.

- We will denote p-value calls made from:
    - ♂ Male only data
    - ♀ Female only data
    - ⚤ All data combined

- We will denote `mp_chooser.json` terms as they are:
    - MALE
    - FEMALE
    - UNSPECIFIED

In the `mp_chooser.json` file each MP term can have different levels:
- Ontology term levels: ABNORMAL, INCREASED, DECREASED.
- Sex levels: FEMALE, MALE, UNSPECIFIED.

*Note:* Sex-specific MP terms, e.g. those with FEMALE and MALE sex levels, and UNSPECIFIED, are never encountered together for the same parameter. In other words, it is either a FEMALE/MALE term or an UNSPECIFIED term available in the `mp_chooser.json` file.

MP term assignment logic can be seen below:

```mermaid
%%{
    init: {
        "theme": "default",
        "themeVariables": {
            "fontSize": "15px"
        },
        "sequence": {
            "useMaxWidth": false
        }
    }
}%%
graph TD;
    Start{Which method is used for the analysis?} --> |MM| MM[Prioritise INCREASED/DECREASED MP term, otherwise use ABNORMAL] --> A
    Start --> |FE or RR| FE_RR[Only use ABNORMAL MP term] --> A

    A{"Is sex-specific MP term available in the mp_chooser file?"}
    A --> |Yes| A2[Use FEMALE/MALE term] --> B{Is ♀ or ♂ call observed?}
    A --> |No| A3[Use UNSPECIFIED term] --> B

    B --> |Yes| B1[Drop ⚤ call and report MP term for ♀ or ♂ call]
    B --> |No| B2[Report ⚤ call]
