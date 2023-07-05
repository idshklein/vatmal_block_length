pacman::p_load(tidyverse,sf,sfnetworks,tidygraph,lwgeom,mapview,nngeo,qgisprocess,cluster,furrr,dbscan,ggrepel,furrr)
qgis_configure()
df <- st_read("Tml_BuiltCells_160223_v3.gpkg")
pols <- df %>%  
  group_by(pl_id) %>% 
  summarise() %>% 
  st_cast("MULTIPOLYGON") %>% 
  st_remove_holes() %>% 
  st_cast("POLYGON") 
ch <- pols %>% 
  mutate(rn = row_number()) %>%{
    q <- (.)
    qgis_run_algorithm("qgis:minimumboundinggeometry",
                       INPUT = q,
                       FIELD = "rn",
                       TYPE = "Convex Hull")
  } %>% `$`("OUTPUT") %>% st_read()
mor <- pols %>% 
  mutate(rn = row_number()) %>%{
    q <- (.)
    qgis_run_algorithm("qgis:minimumboundinggeometry",
                       INPUT = q,
                       FIELD = "rn",
                       TYPE = "Minimum Oriented Rectangle")
  } %>% `$`("OUTPUT") %>% st_read()
ch_ar <- st_area(ch) %>% as.numeric()
mor_ar <- st_area(mor) %>% as.numeric()
pols_ar <- st_area(pols) %>% as.numeric()
tibble(ch_ar,mor_ar,pols_ar) %>% 
  ggplot(aes(x=pols_ar,y=ch_ar)) + 
  geom_point()
tibble(ch_ar,mor_ar,pols_ar) %>% 
  ggplot(aes(x=pols_ar,y=mor_ar)) + 
  geom_point()
tibble(ch_ar,mor_ar,pols_ar) %>% 
  ggplot(aes(x=pols_ar/ch_ar)) + 
  stat_ecdf()
tibble(ch_ar,mor_ar,pols_ar) %>% 
  ggplot(aes(x=pols_ar/mor_ar)) + 
  stat_ecdf()
pols[which(pols_ar/ch_ar > 0.8& pols_ar/ch_ar < 0.9),][1:25,] %>% rowid_to_column("rn") %>% group_split(rn) %>% map(~.x %>% ggplot() + geom_sf())
calculate_segment_angles_directions <- function(geometry) {
  if (st_is(geometry, "LINESTRING")) {
    # Extract the vertices from the linestring
    vertices <- st_coordinates(geometry)
    
    # Calculate the differences between consecutive vertices
    diffs <- diff(vertices)
    
    # Calculate the angles between consecutive segments using atan2
    angles <- atan2(diffs[c(nrow(diffs),1:(nrow(diffs)-1)), "Y"], diffs[c(nrow(diffs),1:(nrow(diffs)-1)), "X"]) - atan2(diffs[, "Y"], diffs[, "X"])
    angles <- angles * 180 / pi
    
    # Adjust the angles to be between 0 and 360 degrees
    angles <- ifelse(angles < 0, angles + 360, angles)
    
    # Calculate the lengths of segments before and after each angle
    lengths_after <- sqrt(diffs[, "X"]^2 + diffs[, "Y"]^2)
    lengths_before <- sqrt(diffs[c(nrow(diffs),1:(nrow(diffs)-1)), "X"]^2 + diffs[c(nrow(diffs),1:(nrow(diffs)-1)), "Y"]^2)
    
    # Calculate the directions based on the angles
    directions <- ifelse(angles > 180, "Left","Right")
    angles <- ifelse(angles > 180, 360 - angles, angles)
    
    # Return a data frame with the angles, lengths before, lengths after, and directions
    return(data.frame(angles = angles, 
                      lengths_before = lengths_before, 
                      lengths_after = lengths_after, 
                      x= vertices[-nrow(vertices),1],
                      y= vertices[-nrow(vertices),2],
                      directions = directions))
  }
  
  if (st_is(geometry, "POLYGON")) {
    # Extract the exterior ring of the polygon
    ring <- st_cast(geometry, "LINESTRING")
    
    # Calculate the angles, lengths before, lengths after, and directions using the linestring calculation
    result <- calculate_segment_angles_directions(ring)
    
    # Return the result
    return(result)
  }
  
  # Return NULL if the geometry type is not supported
  return(NULL)
}

vect <- which(pols_ar > 2000)
pols2 <- pols[vect,]
ch2 <- ch[vect,]
plan(multisession, workers = 8)
res <- future_map(1:nrow(pols2),~{
  line <-pols2[.x,]
  result <- calculate_segment_angles_directions(line) 
  result1 <- result%>%
    mutate(cumlen = cumsum(lengths_before),
           changed_dir = directions != lag(directions),
           changed_dir = ifelse(is.na(changed_dir),F,changed_dir),
           leg = cumsum(changed_dir) +1) %>%
    group_by(leg) %>%
    mutate(len = sum(lengths_before),
           sum_angle = sum(angles),
           directions = ifelse(len > 26&sum_angle >40,directions,NA)) %>%
    ungroup() %>%
    fill(directions,.direction="updown") %>%
    mutate(changed_dir = directions != lag(directions),
           changed_dir = ifelse(is.na(changed_dir),F,changed_dir),
           leg = cumsum(changed_dir) +1) %>%
    group_by(leg) %>%
    mutate(len = sum(lengths_before),
           sum_angle = sum(angles)) %>%
    rowid_to_column("rn") 
  dbs <- dbscan(ungroup(result1[,c("x","y")]),5,1)$cluster
  result1 %>% 
    ungroup() %>%
    mutate(dbs=dbs) %>% 
    st_as_sf(coords = c("x","y"),crs = 2039) %>%
    {
      eds <- data.frame(from = (.)$rn, to = (.)$rn[c(2:nrow(.),1)])
      tbl_graph((.),eds) %>% 
        as_sfnetwork() %>% 
        convert(to_spatial_explicit) %E>% 
        mutate(from_dir = .N()$directions[from],
               to_dir = .N()$directions[to]) %>% 
        filter(from_dir == to_dir) %N>% 
        mutate(comp = group_components()) %>% 
        as_tibble() 
    } %>%
    select(-.tidygraph_node_index) %>% 
    st_drop_geometry() %>% 
    group_by(comp) %>% 
    mutate(comp_sum_len = sum(lengths_after ))%>% 
    ungroup() %>% 
    mutate(samedbscum = cumsum(dbs != lag(dbs,default = last(dbs)))==0) %>% 
    group_split(samedbscum)  %>% 
    bind_rows() %>%
    group_by(dbs) %>%
    mutate(sum_angle2 = sum(angles),
           sum_len2 = sum(lengths_after[1:max(n()-1,1)]),
           cumsum_len2 = cumsum(lengths_after),
           ratio_len2 = cumsum_len2/sum_len2,
           delta = abs(0.5 - ratio_len2),
           slct = delta == min(delta)) %>%
    ungroup() %>%
    pull(comp) %>% 
    unique()
  # mutate(comp2 = ifelse(lag(comp) != comp, lag(comp),ifelse(lead(comp) != comp,lead(comp),comp))) %>% 
  # group_by(slct) %>% 
  # mutate(slct2 = rank(-sum_angle2)) %>% 
  # as.data.frame() %>%
  # filter(slct) %>%
  # arrange(slct2) %>% 
  # head(55) %>%
  # print()
  # result1 %>%
  #   ggplot(aes(x,y,label = rn,color = as.factor(dbs))) +
  #   geom_label_repel(max.overlaps= 50) + 
  #   geom_point()
})
res %>% map_dbl(~.x %>% length) %>% table()
pols2[160,] %>% plot()
reg_anal <- res %>% map_dbl(~.x %>% length) %>% `<=`(2) %>% which()
reg_anal2 <-  res %>% map_dbl(~.x %>% length) %>% `==`(2) %>% which()
ch_anal <- (res %>% map_dbl(~.x %>% length) %>% `>`(2) %>% which())[(ch_ar/pols_ar)[vect][res %>% map_dbl(~.x %>% length) %>% `>`(2)  %>% which()]<1.5] 
pols3 <- pols2[reg_anal,]
ch3 <- ch2[ch_anal,]
plan(multisession, workers = 8)
res2 <- future_map(1:nrow(pols3),~{
  line <-pols3[.x,]
  result <- calculate_segment_angles_directions(line) 
  result1 <- result%>%
    mutate(cumlen = cumsum(lengths_before),
           changed_dir = directions != lag(directions),
           changed_dir = ifelse(is.na(changed_dir),F,changed_dir),
           leg = cumsum(changed_dir) +1) %>%
    group_by(leg) %>%
    mutate(len = sum(lengths_before),
           sum_angle = sum(angles),
           directions = ifelse(len > 26&sum_angle >40,directions,NA)) %>%
    ungroup() %>%
    fill(directions,.direction="updown") %>%
    mutate(changed_dir = directions != lag(directions),
           changed_dir = ifelse(is.na(changed_dir),F,changed_dir),
           leg = cumsum(changed_dir) +1) %>%
    group_by(leg) %>%
    mutate(len = sum(lengths_before),
           sum_angle = sum(angles)) %>%
    rowid_to_column("rn") 
  dbs <- dbscan(ungroup(result1[,c("x","y")]),5,1)$cluster
  result1 %>% 
    ungroup() %>%
    mutate(dbs=dbs) %>% 
    st_as_sf(coords = c("x","y"),crs = 2039) %>%
    {
      eds <- data.frame(from = (.)$rn, to = (.)$rn[c(2:nrow(.),1)])
      tbl_graph((.),eds) %>% 
        as_sfnetwork() %>% 
        convert(to_spatial_explicit) %E>% 
        mutate(from_dir = .N()$directions[from],
               to_dir = .N()$directions[to]) %>% 
        filter(from_dir == to_dir) %N>% 
        mutate(comp = group_components()) %>% 
        as_tibble() 
    } %>%
    select(-.tidygraph_node_index) %>% 
    st_drop_geometry() %>% 
    group_by(comp) %>% 
    mutate(comp_sum_len = sum(lengths_after ))%>% 
    ungroup() %>% 
    mutate(samedbscum = cumsum(dbs != lag(dbs,default = last(dbs)))==0) %>% 
    group_split(samedbscum)  %>% 
    bind_rows() %>%
    group_by(dbs) %>%
    mutate(sum_angle2 = sum(angles),
           sum_len2 = sum(lengths_after[1:max(n()-1,1)]),
           cumsum_len2 = cumsum(lengths_after),
           ratio_len2 = cumsum_len2/sum_len2,
           delta = abs(0.5 - ratio_len2),
           slct = delta == min(delta),
           ns=n()) %>%
    ungroup() %>%
    mutate(complag2 = lag(comp,default = last(comp))==2,
           complead2 = lead(comp,default = first(comp))==2,
           comp2=ifelse(complag2|complead2,2,comp)) %>%
    group_by(slct) %>%
    mutate(slct2 = rank(-sum_angle2)) %>%
    as.data.frame() %>%
    arrange(slct2)
})
res2 %>% 
  imap(~.x %>% mutate(idx = .y)) %>% 
  bind_rows() %>% 
  group_by(idx) %>% 
  mutate(comt = max(comp2)) %>% 
  slice(1) %>% 
  filter(comt==2, round(sum_angle2) == 180,ns==2) %>% 
  pull(idx) %>% 
  `[`(93)

pols3[904,] %>% st_write("wert.gpkg",delete_layer = T)
res2[[904]] %>%
  mutate(morethan50 = sum_angle2 > 50) %>% 
  arrange(rn) %>% 
  mutate(slct3= slct&morethan50,
         cumm = cumsum(slct3)) %>% 
  select(cumm, everything()) %>% 
  group_split(cumm) %>% 
  {
    q <- (.)
    if(q[[1]]$cumm[1] == 0){
      lt <- length(q)
      ret <- q[c(2:lt,1)]  
    }else{
      ret <- q
    }
    q
  } %>% 
  bind_rows() %>% 
  mutate(cumm = cummax(cumm)  )%>%
  group_by(cumm) %>% 
  summarise(lengths_after=sum(lengths_after)) 
pols3[1000,]%>%st_cast("POINT") %>%  mapview()
res2[res2 %>% 
  imap(~.x %>% mutate(idx = .y)) %>% 
  bind_rows() %>% 
  group_by(idx) %>% 
  mutate(comt = max(comp2)) %>% 
  slice(1) %>% 
  filter(comt==1, round(sum_angle2) == 180) %>%
  pull(idx)] %>% 
  imap(~.x %>% mutate(idx = .y)) %>% 
  bind_rows() %>% 
  group_by(idx) %>% 
  slice(1) %>% 
  pull(ns) %>%
  `==`(366) %>% 
  which()
res2[[3684]] 
pols3[167,] %>%st_cast("POINT") %>% rowid_to_column("rn") %>% 
    ggplot(aes(x,y,label = rn,color = as.factor(dbs))) +
    geom_label_repel(max.overlaps= 50) +
    geom_point()
  ggplot(aes(x=sum_len2,y=sum_angle2)) + 
  geom_point() + 
  coord_cartesian(xlim = c(0,50),ylim = c(0,800))
  ggplot(aes(x=factor(comt) ,y =sum_angle2)) + 
  geom_violin()
res2[[1]] %>% as.data.frame() %>% head(100)
res2 %>% as_tibble() %>% 
  filter(value<500) %>% 
  ggplot(aes(value)) + 
  stat_ecdf()
