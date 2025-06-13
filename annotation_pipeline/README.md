# IMPC Annotation Pipeline
The IMPC annotation pipeline assigns Mammalian Phenotype (MP) terms to significant genetic effects based on a p-value threshold of 0.0001. The goal is to associate phenotypic observations with the corresponding genetic modifications.
At the IMPC, genetic effects are identified using three statistical analysis platforms:
1. Linear Mixed Model framework or MM
2. Fisher Exact Test framework or FE
3. Reference Range Plus Test framework or RR

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
    A[Repeat separately for Overall, Females and Males] --> B{Is genotype effect significant?}
    B -- No --> C[Do not assign MP term]
    B -- Yes --> D{Is the direction of genotype effect specified?}
    D -- No --> E[Select Abnormal term]
    D -- Yes --> F{Is there any conflict of direction? <br> example: Low.Decrease/High.Increase <br> or Male.Decrease/Female.Increase}
    F -- No --> G[Choose an MP term corresponding to the direction of change]
    F -- Yes --> H[Select Abnormal term]
```