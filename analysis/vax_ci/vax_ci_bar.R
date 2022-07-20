library(tidyverse)
library(RColorBrewer)
library(htmltools)


fs::dir_create(here::here("output", "vax_ci"))



data <- readr::read_delim(
  "analysis/vax_ci/opensafely20211208", 
  delim = "\t",
  escape_double = FALSE,
  col_names = FALSE,
  skip = 1
  ) %>%
  rename(group=X1) %>%
  # rename(group=X1,dose1=X2,dose2=X3,dose3=X4) %>%
  mutate(across(starts_with("X"), 
                ~as.numeric(str_extract(.x, "\\d+\\.\\d+"))/100)) %>%
  transmute(
    group,
    "unvaccinated" = 1-X2,
    "1 dose" = X2 - X3,
    "2 doses" = X3 - X4,
    "3 doses" = X4
    ) %>%
  filter(!(group %in% c("16-17", "care home", "LD (aged 16-64)"))) %>%
  mutate(across(group, 
                ~if_else(
                  str_detect(.x, "shielding"),
                  "16-69 years\n& shielding",
                  str_c(.x, " years")
                )))

groups <- c("unvaccinated", "1 dose", "2 doses", "3 doses")

pal <- c("grey", brewer.pal(n=3, name="Dark2"))
names(pal) <- groups

data %>%
  mutate(across(group, factor, levels = data$group)) %>%
  pivot_longer(cols = -group) %>%
  mutate(across(name, factor, levels = groups)) %>%
  ggplot(aes(x = reorder(group,-as.integer(group)), y = value, fill = name)) +
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE)
    ) +
  scale_fill_manual(
    name = NULL,
    values = pal
    ) +
  scale_x_discrete(
    name = NULL
  ) +
  scale_y_continuous(
    name = NULL,
    labels = scales::percent_format(accuracy = 1),
    expand = c(0.01,0.01)
    ) +
  labs(
    title = "Priority groups split by vaccine status as of 8 December 2021",
    caption = "Data avaiable from the OpenSAFELY vaccine coverage report archive: https://tinyurl.com/os-vax-211208"
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.caption.position = "plot",
    plot.caption = element_text(size=8),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(l=0,t=0,r=0,b=15))
  )

ggsave(
  filename = here::here("output", "vax_ci", "vax_ci_bar.png"),
  width = 14, height = 14, units = "cm"
)
