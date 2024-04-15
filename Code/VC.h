#pragma once
#include "Lexer.h"

// GLOBAL
extern struct DFA *vc_lexer_automata;

void               vc_lexer_automata_init();
static inline void vc_lexer_automata_free() { dfa_free(vc_lexer_automata); }

