library(DiagrammeR)

internal_lab <- c(
  "State-dependent?",
  "Interruptable?",
  "Interruptable?",
  "Interruption \nstate-dependent?",
  "Interruption \nstate-dependent?"
)

terminal_lab <- c(
  "6",
  "5",
  "4",
  "3",
  "2",
  "1"
)

grViz("
  digraph delay_types {
    
    node [shape = box, style = filled, fillcolor = lightsteelblue1]
    a [label = '@@1-1']
    b [label = '@@1-2']
    c [label = '@@1-3']
    d [label = '@@1-4']
    e [label = '@@1-5']
    
    node [shape = circle, fixedsize = true, style = filled, fillcolor = LightGray]
    f [label = '@@2-1']
    g [label = '@@2-2']
    h [label = '@@2-3']
    i [label = '@@2-4']
    j [label = '@@2-5']
    k [label = '@@2-6']
    
    a->b [label = 'Yes']
    a->c [label = 'No']
    b->d [label = 'Yes']
    b->h [label = 'No']
    d->f [label = 'Yes']
    d->g [label = 'No']
    c->e [label = 'Yes']
    c->k [label = 'No']
    e->i [label = 'Yes']
    e->j [label = 'No']
    
  }
  
  [1]: internal_lab
  [2]: terminal_lab
")