#pragma once
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define STRUCT_SIZE(sname, mname, nelem) MAX(sizeof(sname), offsetof(sname, mname) + (nelem) * sizeof(*((sname){}).mname))

#define ASCII_TABLE_SIZE 128

typedef uint64_t token;
#define TOKEN_IDENTIFIER_BIT (((token)1)<<0)
#define TOKEN_KEYWORD_BIT    (((token)1)<<1)
#define TOKEN_OPERATOR_BIT   (((token)1)<<2)
#define TOKEN_SEPERATOR_BIT  (((token)1)<<3)
#define TOKEN_LITERAL_BIT    (((token)1)<<4)

// Things that are ignored by the lexer
#define TOKEN_IGNORE_BIT     (((token)1)<<5)

static const char *token_names[] = {
    "IDENTIFIER",
    "KEYWORD",
    "OPERATOR",
    "SEPERATOR",
    "LITERAL"
};
#define NUM_TOKEN_KIND (sizeof(token_names) / sizeof(*token_names))
static_assert(NUM_TOKEN_KIND < 64, "too many token kinds");


// Used to mark a state as terminal without defining the token(s) for it
#define TOKEN_UNDEFINED      (((token)1)<<63)
// Used to mark a state as nonterminal
#define NONTERMINAL          ((token)0)

#define EPSILON 0

#define DIGIT  "0123456789"
#define LETTER "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_"

typedef char input;

struct Transition {
    input character;
    size_t successor;
};
static inline bool transition_eq(const struct Transition *lhs, const struct Transition *rhs) { return lhs->successor == rhs->successor && lhs->character == rhs->character; }

// TODO: Rename this to NFAState
struct State {
    token terminality;
    size_t num_transition;
    struct Transition transitions[];
};

// Initialize a state
struct State *state_init                (token terminality, size_t num_transition);
// Like state_init but specify the transitions right aways
struct State *state_create              (token terminality, size_t num_transition, input in[num_transition], size_t succs[num_transition]);
// Copy a state, add an `extra_transitions` number of transitions to the start of the transtions array that can be writen to
// old state:
// | N transitions |
// new state:
// | extra transitions (undefined) | N old transitions |
struct State *state_copy                (const struct State *state, size_t extra_transitions);
// Resize a state (truncate the last `num_deficit` transtions)
struct State *state_mut_shrink          (struct State *state, size_t num_deficit);
// Sort the transtions based on input character and successor state
void          state_mut_sort_transition (struct State *state);
// Remove transtions that except character
void          state_mut_strip_transition(struct State *state, input character);
// Strip duplicated transitions, invoke state_mut_sort_transition to efficiently detect duplication
void          state_mut_strip_duplicate (struct State *state);

static inline void state_free(struct State *state) { free(state); }

struct NFA {
    // The number of state
    size_t num_state;
    // The number of starting states
    size_t num_initial;
    // An array of states
    // NOTE: The first num_initial states are the starting state
    struct State *states[];
};

// Some test to check if automata is a valid NFA 
static inline void
nfa_check_correctness(const struct NFA *automata) {
    assert(automata->num_initial <= automata->num_state);
    for (size_t i = 0; i < automata->num_state; ++i)
        for (size_t t = 0; t < automata->states[i]->num_transition; ++t) {
            assert(automata->states[i]->transitions[t].successor < SIZE_MAX);
            assert(automata->states[i]->transitions[t].successor < automata->num_state);
    }
}

// Initialize NFA with num_state states in which the first num_initial states are starting states
struct NFA *nfa_init                        (size_t num_state, size_t num_initial);
void        nfa_free                        (struct NFA *automata);
// Set all TOKEN_UNDEFINED to definition. Only affect states where terminality == TOKEN_UNDEFINED
void        nfa_define_token                (struct NFA *automata, token definition);
// print representation to stream
void        nfa_fprint                      (FILE *stream, const struct NFA *automata);
// Common constructions
// nfa_verbatim will match the `pattern` string exactly
struct NFA *nfa_verbatim                    (const input *pattern);
// nfa_any will match any characters in the `options` string
struct NFA *nfa_any                         (const input *options);
// union of two NFAs
struct NFA *nfa_union                       (const struct NFA *lhs, const struct NFA *rhs);
struct NFA *nfa_union_many                  (const struct NFA *automata, ...);
struct NFA *nfa_optional                    (const struct NFA *automata);
// kleene start of NFA
struct NFA *nfa_kleene                      (const struct NFA *automata);
// concat of two NFAs
struct NFA *nfa_concat                      (const struct NFA *prefix, const struct NFA *suffix);
struct NFA *nfa_concat_many                 (const struct NFA *first, ...);
// Attempt NFA minimalization
void        nfa_mut_minimize                (struct NFA *automata);

// Maybe remove these
void        nfa_mut_eliminate_epsilon       (struct NFA *automata);
void        nfa_mut_strip_unreachable_states(struct NFA *automata);
void        nfa_mut_strip_dead_states       (struct NFA *automata);
void        nfa_mut_remove_duplicate_states (struct NFA *automata);

// Common NFA construction
static inline struct NFA *
nfa_digit() {
    return nfa_any(DIGIT);
}

static inline struct NFA *
nfa_letter() {
    return nfa_any(LETTER);
}

static inline struct NFA *
nfa_alphanum() {
    return nfa_any(DIGIT LETTER);
}

struct DFAState {
    token terminality;
    size_t num_transitions;
    struct Transition unique_transitions[]; // NOTE: the entries are ordered so binary search can be performed
};

struct DFAState    *dfa_state_init(token terminality, size_t num_transitions);
static inline void  dfa_state_free(struct DFAState *state) { free(state); }
size_t              dfa_state_search(const struct DFAState *state, input in);

struct DFA {
    size_t num_states;
    struct DFAState *states[];    
};

// Construct a deterministic finite state automata from a nondeterministic finite state automata
// The NFA must be minimized before begin passed to this function
struct DFA *dfa_powerset_construction(const struct NFA *nfa);
// Release resource held by dfa
void        dfa_free                 (struct DFA *dfa);
// Print a representation of a dfa into stream
void        dfa_fprint               (FILE *stream, const struct DFA *automata);
// Feed a NULL terminated string to dfa and return the resulting token (NONTERMINAL if is not recoginized as a token)
token       dfa_recognize            (const struct DFA *automata, input *str);
// Feed a NFA one character at a time (NOTE: The starting state is 0, and there exist a dead state with value SIZE_MAX)
size_t      dfa_consume              (const struct DFA *automata, size_t state, input in);
// Retreive the token of this state (NONTERMINAL if is non terminal state)
token       dfa_state_token          (const struct DFA *automata, size_t state);
// TODO: Generate C code
