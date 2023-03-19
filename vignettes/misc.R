## Saturation curve 

Saturation curve of average beta-diversity between communities with
community 1.

In this part, we demonstrate that the average distance from
other communities to community 1 will reach to a threshold of
nutrients, after which the average distance won't increase along with
the total concentration of nutrients.

Let us first define a function calculating the mean distance to the
first community.



## UMAP distance saturation analysis

This shows how distance saturation could be calculated. Not evaluated
currently.

```{r df, eval=TRUE}
distance_saturation_data <- data.frame(concentration = integer(),
                                       community = integer(),
                                       average_distance = numeric())

#umap_CRM_coor <- cbind(umap_CRM_coor, concentration, medium, community)
umap_CRM_coor <- data.frame(reducedDim(tse), colData(tse)) # cbind(umap_CRM_coor, concentration, medium, community)
for (res_conc_type in unique(concentration)){
    for (com_type in unique(community)){
        ave_dist <- average_distance(umap_CRM_coor, res_conc_type, com_type)
        distance_saturation_data[nrow(distance_saturation_data)+1,] <-
            c(res_conc_type, com_type, ave_dist)
    }
}
# View(distance_saturation_data)
distance_saturation_data$average_distance <- as.numeric(distance_saturation_data$average_distance)
distance_saturation_data$concentration <- as.factor(distance_saturation_data$concentration)
distance_saturation_data$community <- as.factor(distance_saturation_data$community)

# distance_saturation_data_plot
p <- ggplot(distance_saturation_data,
                                 aes(concentration, average_distance,
                                     color = community,
                                     group = community)) +
    geom_line() +
    geom_point() +
    scale_shape_manual(values = c(0, 1, 2, 5, 6, 8, 15, 16, 17, 18)) +
    labs(x = "Resource concentration",
         y = "Average distance between communities in UMAP") +
    theme_bw()

print(p)
```
