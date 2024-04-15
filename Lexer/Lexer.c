#include "Lexer.h"
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>



struct State *
state_init(token terminality, size_t num_transition) {
    struct State *state = malloc(STRUCT_SIZE(struct State, transitions, num_transition));
    *state = (struct State){.terminality=terminality,.num_transition=num_transition};
    return state;
}

struct State *
state_create(token terminality, size_t num_transition, input in[num_transition], size_t succs[num_transition]) {
    struct State *state = state_init(terminality, num_transition);
    for (size_t i = 0; i < num_transition; ++i) 
        state->transitions[i] = (struct Transition){.character=in[i],.successor=succs[i]};
    return state;
}

struct State *
state_copy(const struct State *state, size_t extra_transitions) {
    struct State *copy = state_init(state->terminality, state->num_transition + extra_transitions);
    memcpy(copy->transitions + extra_transitions, state->transitions, state->num_transition * sizeof(struct Transition));
    return copy;
}

struct State *
state_mut_shrink(struct State *state, size_t num_deficit) {
    assert(state->num_transition >= num_deficit);
    state->num_transition -= num_deficit;
    struct State *new = realloc(state, STRUCT_SIZE(struct State, transitions, state->num_transition));
    return new;
}

void 
state_mut_sort_transition(struct State *state) {
    // NOTE: Implement a more efficient algo
    bool sorted = false;
    while (!sorted) {
        sorted = true;
        for (size_t i = 1; i < state->num_transition; ++i) {
            if (state->transitions[i - 1].character < state->transitions[i].character) continue;
            if (state->transitions[i - 1].character == state->transitions[i].character &&
                    state->transitions[i - 1].successor <= state->transitions[i].successor)
                continue;

            sorted = false;
            struct Transition temp = state->transitions[i - 1];
            state->transitions[i - 1] = state->transitions[i];
            state->transitions[i] = temp;
        }
    }
}

void 
state_mut_strip_transition(struct State *state, input character) {
    struct Transition *read, *write;
    for (read = state->transitions, write = state->transitions;
            read < state->transitions + state->num_transition;
            ++read) {
        if (read->character != character)
            *(write++) = *read;
    }
    state->num_transition = write - state->transitions;
}

void
state_mut_strip_duplicate(struct State *state) {
    if (state->num_transition == 0) return;
    state_mut_sort_transition(state);
    struct Transition *read, *write;
    for (read = state->transitions + 1, write = read;
            read < state->transitions + state->num_transition;
            ++read) {
        if (!transition_eq(read, write - 1)) 
            *(write++) = *read;
    }
    state->num_transition = write - state->transitions;
}



struct NFA *
nfa_init(size_t num_state, size_t num_initial) {
    assert(num_initial <= num_state);
    struct NFA *automata = malloc(STRUCT_SIZE(struct NFA, states, num_state));
    *automata = (struct NFA){.num_state = num_state, .num_initial = num_initial};
    memset(automata->states, 0, sizeof(struct State *) * automata->num_state);
    return automata;
}

void
nfa_free(struct NFA *automata) {
    for (size_t i = 0; i < automata->num_state; ++i)
        state_free(automata->states[i]);
    free(automata);
}

void 
nfa_define_token(struct NFA *automata, token definition) {
    for (size_t i = 0; i < automata->num_state; ++i) {
        struct State *state = automata->states[i];
        if (state->terminality == TOKEN_UNDEFINED)
            state->terminality = definition;
    }
}

void
nfa_fprint(FILE *stream, const struct NFA *automata) {
    struct State *const *state = automata->states;
    for (size_t idx = 0; idx < automata->num_state; ++idx, ++state) {
        fprintf(stream, "==[STATE %lu]\n", idx);
        if (idx < automata->num_initial)
            fprintf(stream, "  STARTING\n");
        if ((*state)->terminality != NONTERMINAL) {
            fputs("  TERMINAL: ", stream);
            if ((*state)->terminality & TOKEN_OPERATOR_BIT)
                fputs("OPERATOR ", stream);
            if ((*state)->terminality & TOKEN_LITERAL_BIT)
                fputs("LITERAL ", stream);
            if ((*state)->terminality & TOKEN_IDENTIFIER_BIT)
                fputs("IDENTIFIER ", stream);
            if ((*state)->terminality & TOKEN_KEYWORD_BIT)
                fputs("KEYWORD ", stream);
            if ((*state)->terminality & TOKEN_SEPERATOR_BIT)
                fputs("SEPERATOR ", stream);
            fputc('\n', stream);
        }
        for (size_t t = 0; t < (*state)->num_transition; ++t) {
            struct Transition *transition = (*state)->transitions + t;
            if (isprint(transition->character)) {
                fprintf(stream, "  '%c' -> %lu\n", transition->character, transition->successor);
            } else if (transition->character == 0) {
                fprintf(stream, "  epsilon -> %lu\n", transition->successor);
            } else {
                fprintf(stream, "  '\\%x' -> %lu\n", transition->character, transition->successor);
            }
        }
    }
}

// Input:
// - pattern : null terminated string
// Output:
// Nondeterministic automata
struct NFA *
nfa_verbatim(const input *pattern) {
    size_t vertex_num = strlen(pattern);
    struct NFA *automata = nfa_init(1 + vertex_num, 1);
    for (size_t i = 0; i < vertex_num; ++i) {
        assert(isascii(pattern[i]));
        automata->states[i] = state_create(NONTERMINAL, 1, (input[]){pattern[i]}, (size_t[]){i + 1});
    }
    automata->states[vertex_num] = state_create(TOKEN_UNDEFINED, 0, NULL, NULL);

    nfa_check_correctness(automata);
    return automata;
}

struct NFA *
nfa_any(const input *options) {
    size_t num_options = strlen(options);
    struct NFA *automata = nfa_init(2, 1);
    
    automata->states[0] = state_init(NONTERMINAL, num_options);
    for (size_t i = 0; i < num_options; ++i)
        automata->states[0]->transitions[i] = (struct Transition){.character=options[i], .successor=1};
    automata->states[1] = state_init(TOKEN_UNDEFINED, 0);

    return automata;
}

struct NFA *
nfa_union(const struct NFA *lhs, const struct NFA *rhs) { return nfa_union_many(lhs, rhs, NULL); }

struct NFA *
nfa_union_many(const struct NFA *first, ...) {
    va_list ap;
    size_t num_state   = 0,
           num_initial = 0;

    va_start(ap, first);
    for (const struct NFA *other = first; other != NULL; other = va_arg(ap, const struct NFA *)) {
        num_state   += other->num_state;
        num_initial += other->num_initial;
    }
    va_end(ap);

    struct NFA *result = nfa_init(num_state, num_initial);
    size_t initial_state_offset    = 0,
           noninitial_state_offset = num_initial;
    va_start(ap, first);
    for (const struct NFA *other = first; other != NULL; other = va_arg(ap, const struct NFA *)) {
        assert(initial_state_offset < num_initial);
        assert(noninitial_state_offset < num_state);
        for (size_t i = 0; i < other->num_initial; ++i) {
            struct State *state = state_copy(other->states[i], 0);
            for (struct Transition *trans = state->transitions;
                    trans < state->transitions + state->num_transition; 
                    ++trans) {
                trans->successor = trans->successor < other->num_initial ? trans->successor + initial_state_offset : trans->successor + noninitial_state_offset - other->num_initial;
            }
            assert(initial_state_offset + i < num_initial);
            result->states[initial_state_offset + i] = state;
        }
        for (size_t i = other->num_initial; i < other->num_state; ++i) {
            struct State *state = state_copy(other->states[i], 0);
            for (struct Transition *trans = state->transitions;
                    trans < state->transitions + state->num_transition; 
                    ++trans) {
                trans->successor = trans->successor < other->num_initial ? trans->successor + initial_state_offset : trans->successor + noninitial_state_offset - other->num_initial;
            }
            assert(noninitial_state_offset + i - other->num_initial < num_state);
            result->states[noninitial_state_offset + i - other->num_initial] = state;
        }
        initial_state_offset += other->num_initial;
        noninitial_state_offset += other->num_state - other->num_initial;
    }
    assert(initial_state_offset == num_initial);
    assert(noninitial_state_offset == num_state);
    va_end(ap);

    return result;
}

struct NFA *
nfa_kleene(const struct NFA *automata) {
    size_t num_state = automata->num_state,
           num_initial = automata->num_initial;
    struct NFA *kleene = nfa_init(num_state, num_initial);
    for (size_t i = 0; i < num_state; ++i) {
        struct State *state;
        if (automata->states[i]->terminality != NONTERMINAL) {
            state = state_copy(automata->states[i], num_initial);
            for (size_t j = 0; j < num_initial; ++j)
                state->transitions[j] = (struct Transition){.character=EPSILON, .successor=j};
        } else {
            state = state_copy(automata->states[i], 0);
        }
        if (i < num_initial)
            state->terminality = TOKEN_UNDEFINED;
        kleene->states[i] = state;
    }
    nfa_check_correctness(kleene);
    return kleene;
}

struct NFA *
nfa_concat(const struct NFA *prefix, const struct NFA *suffix) {
    return nfa_concat_many(prefix, suffix, NULL);
}


struct NFA *
nfa_concat_many(const struct NFA *first, ...) {
    va_list ap;
    if (first == NULL)
        return nfa_verbatim("");

    size_t num_initial = first->num_initial,
           num_state   = 0;

    va_start(ap, first);
    for (const struct NFA *nfa = first; nfa != NULL; nfa = va_arg(ap, const struct NFA *))
        num_state += nfa->num_state;
    va_end(ap);

    struct NFA *automata = nfa_init(num_state, num_initial);

    va_start(ap, first);
    size_t offset = 0;
    const struct NFA *prefix = first, *suffix = va_arg(ap, const struct NFA *);
    for (; suffix != NULL; prefix = suffix, suffix = va_arg(ap, const struct NFA *)) {
        for (size_t i = 0; i < prefix->num_state; ++i) {
            struct State *state;
            if (prefix->states[i]->terminality == NONTERMINAL) {
                state = state_copy(prefix->states[i], 0);
                for (struct Transition *trans = state->transitions;
                        trans < state->transitions + state->num_transition; 
                        ++trans)
                    trans->successor += offset;
            } else {
                state = state_copy(prefix->states[i], suffix->num_initial);
                for (size_t i = 0; i < suffix->num_initial; ++i) 
                    state->transitions[i] = (struct Transition){.character = EPSILON, .successor = prefix->num_state + i};
                for (struct Transition *trans = state->transitions;
                        trans < state->transitions + state->num_transition; 
                        ++trans)
                    trans->successor += offset;
                state->terminality = NONTERMINAL;
            }
            automata->states[offset + i] = state;
        }
        // updating offset
        offset += prefix->num_state;
    }
    assert(offset + prefix->num_state == num_state);
    for (size_t i = 0; i < prefix->num_state; ++i) {
        struct State *state = state_copy(prefix->states[i], 0);
        for (struct Transition *trans = state->transitions;
                trans < state->transitions + state->num_transition; 
                ++trans)
            trans->successor += offset;
        automata->states[offset + i] = state;
    }

    va_end(ap);

    return automata;
}

struct NFA *
nfa_optional(const struct NFA *automata) {
    struct NFA *new_automata = nfa_init(automata->num_state, automata->num_initial);
    for (size_t i = 0; i < automata->num_state; ++i) {
        new_automata->states[i] = state_copy(automata->states[i], 0);
        if (i < automata->num_initial)
            new_automata->states[i]->terminality = TOKEN_UNDEFINED;
    }
    return new_automata;
}

// No lambda function :(
static size_t
nfa_mut_eliminate_epsilon_fetch_reachable_state(const struct NFA *automata, size_t state, size_t num_reached, bool *visited, size_t *scratch, size_t *non_epsilon, token *terminality) {
    assert(state < automata->num_state);
    assert(num_reached <= automata->num_state);

    if (visited[state]) return num_reached;
    visited[state] = true;
    scratch[num_reached] = state;
    ++num_reached;

    const struct State *s = automata->states[state];
    *terminality |= s->terminality;
    for (const struct Transition *t = s->transitions; t < s->transitions + s->num_transition; ++t) {
        if (t->character != EPSILON) {
            ++(*non_epsilon);
            continue;
        }
        num_reached = nfa_mut_eliminate_epsilon_fetch_reachable_state(automata, t->successor, num_reached, visited, scratch, non_epsilon, terminality);
    }

    return num_reached;
}

void
nfa_mut_eliminate_epsilon(struct NFA *automata) {
    bool   *visited = calloc(automata->num_state, sizeof(bool));
    size_t *scratch = calloc(automata->num_state, sizeof(size_t));
    
    for (size_t idx = 0; idx < automata->num_state; ++idx) {
        struct State *oldstate = automata->states[idx];
        
        visited[idx] = true; // Not counting this state itself
        size_t extra_trans = 0;
        size_t num_reachable = 0;
        token terminality = NONTERMINAL;
        
        for (const struct Transition *it = oldstate->transitions;
                it < oldstate->transitions + oldstate->num_transition; 
                ++it) {
            if (it->character == EPSILON)
                num_reachable = nfa_mut_eliminate_epsilon_fetch_reachable_state(automata, it->successor, num_reachable, visited, scratch, &extra_trans, &terminality);
        }
        // Remove epsilon
        state_mut_strip_transition(oldstate, EPSILON);
        struct State *newstate = state_copy(oldstate, extra_trans);
        newstate->terminality |= terminality;
        state_free(oldstate);

        struct Transition *write = newstate->transitions;
        for (size_t st = 0; st < num_reachable; ++st) {
            struct State *state = automata->states[scratch[st]];
            for (const struct Transition *tr = state->transitions;
                    tr < state->transitions + state->num_transition;
                    ++tr) {
                if (tr->character != EPSILON) {
                    *write = *tr;
                    assert(write->successor < automata->num_state);
                    ++write;
                }
                assert(tr->successor < automata->num_state);
            }
        }
        assert(write - newstate->transitions == extra_trans);
        // assert(newstate->num_transition == newstate->num_transition);


        state_mut_strip_duplicate(newstate);
        for (const struct Transition *tr = newstate->transitions;
                tr < newstate->transitions + newstate->num_transition;
                ++tr)
            assert(tr->successor < automata->num_state);

        automata->states[idx] = newstate;
        memset(visited, 0, sizeof(*visited) * automata->num_state);

    }

    for (size_t i = 0; i < automata->num_state; ++i) {
        const struct State *state = automata->states[i];
        for (const struct Transition *trans = state->transitions; 
                trans < state->transitions + state->num_transition;
                ++trans)
            assert(trans->character != EPSILON);
    }

    free(visited);
    free(scratch);
}

static size_t 
nfa_mut_strip_unreachable_states_dfs(const struct NFA *automata, size_t state, bool *visited, size_t count) {
    assert(state < automata->num_state);
    if (visited[state]) return count;
    visited[state] = true;
    ++count;
    
    const struct State *st = automata->states[state];
    for (const struct Transition *tr = st->transitions; 
            tr < st->transitions + st->num_transition;
            ++tr) {
        count = nfa_mut_strip_unreachable_states_dfs(automata, tr->successor, visited, count);
    }

    return count;

}

void 
nfa_mut_strip_unreachable_states(struct NFA *automata) {
    bool *visited = calloc(automata->num_state, sizeof(bool));
    size_t *remap = calloc(automata->num_state, sizeof(size_t)),
           num_reachable = 0;

    for (size_t i = 0; i < automata->num_initial; ++i)
        num_reachable = nfa_mut_strip_unreachable_states_dfs(automata, i, visited, num_reachable);

    size_t new_name = 0;
    for (size_t old_name = 0; old_name < automata->num_state; ++old_name) {
        if (visited[old_name])
            remap[old_name] = new_name++;
        else 
            remap[old_name] = SIZE_MAX;
    }
    assert(new_name == num_reachable);

    for (size_t state = 0; state < automata->num_state; ++state) {
        if (!visited[state])
            state_free(automata->states[state]);
        else 
            automata->states[remap[state]] = automata->states[state];
    }
    
    for (size_t i = 0; i < num_reachable; ++i) {
        struct State *state = automata->states[i];
        for (struct Transition *trans = state->transitions;
                trans < state->transitions + state->num_transition;
                ++trans) {
            trans->successor = remap[trans->successor];
        }
    }
    automata->num_state = num_reachable;

    free(visited);
    free(remap);
    nfa_check_correctness(automata);
}

void
nfa_mut_strip_dead_states(struct NFA *automata) {
    return; // TODO: Fix
    bool *open = calloc(automata->num_state, sizeof(bool));
    bool stable = false;
    while (!stable) {
        stable = true;
        for (size_t i = 0; i < automata->num_state; ++i) {
            if (open[i])
                continue;

            struct State *state = automata->states[i];
            if (state->terminality != NONTERMINAL) {
                open[i] = true;
                stable = false;
            } else {
                for (struct Transition *trans = state->transitions; 
                        trans < state->transitions + state->num_transition;
                        ++trans) {
                    if (open[trans->successor]) {
                        open[i] = true;
                        stable = false;
                        break;
                    }
                }
            }
        }
    }

    size_t read, write;
    for (read = 0, write = read; 
            read < automata->num_state;
            ++read) {
        if (open[read]) {
            automata->states[write++] = automata->states[read];
        } else {
            state_free(automata->states[read]);
        }
    }
    
    automata->num_state = write;
    size_t num_init = 0;
    for (size_t i = 0; i < automata->num_initial; ++i)
        num_init += open[i] ? 1 : 0;
    automata->num_initial = num_init;
    
    free(open);
    nfa_check_correctness(automata);
}

void
nfa_mut_remove_duplicate_states(struct NFA *automata) {
    // TODO: Implement
}

void
nfa_mut_minimize(struct NFA *automata) {
    nfa_mut_eliminate_epsilon(automata);
    nfa_mut_strip_dead_states(automata);
    nfa_mut_strip_unreachable_states(automata);
    // nfa_mut_remove_duplicate_states(automata);

    nfa_check_correctness(automata);
}


// DFAState functions
struct DFAState *
dfa_state_init(token terminality, size_t num_transitions) {
    struct DFAState *state = malloc(STRUCT_SIZE(struct DFAState, unique_transitions, num_transitions));
    *state = (struct DFAState){terminality, num_transitions};
    return state;
}

size_t
dfa_state_search(const struct DFAState *state, input in) {
    for (size_t i = 1; i < state->num_transitions; ++i)
        assert(state->unique_transitions[i - 1].character < state->unique_transitions[i].character);

    const struct Transition *arr = state->unique_transitions;
    size_t num_trans = state->num_transitions;

    while (num_trans) {
        const struct Transition *middle = arr + num_trans / 2;
        if (middle->character == in) {
            return middle->successor;
        } else if (middle->character < in) {
            arr += num_trans / 2 + 1;
            num_trans -= num_trans / 2 + 1;
        } else { // character > in
            num_trans /= 2;
        }
    }
    return SIZE_MAX;
}


static inline bool 
get_bit(const uint64_t *bitmap, size_t bit) { return (bool)((bitmap[bit / 64] >> (bit % 64)) & 1); }

static inline void
set_bit(uint64_t *bitmap, size_t bit, bool value) {
    if (value) {
        bitmap[bit / 64] |= (uint64_t)1 << (bit % 64);
    } else {
        bitmap[bit / 64] &= !((uint64_t)1 << (bit % 64));
    }
}

static inline bool
bitmap_eq(const uint64_t *lhs, const uint64_t *rhs, size_t nword) {
    for (size_t i = 0; i < nword; ++i) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }
    return true;
}

static inline bool
bitmap_iszero(const uint64_t *bitmap, size_t nword) {
    for (size_t i = 0; i < nword; ++i)
        if (bitmap[i] != 0) return false;
    return true;
}


struct PowersetMap {
    const size_t bitmap_size; // In number of uint64_t elements
    size_t size, capacity;
    uint64_t *keys;
    struct DFAState **states;
};

static struct PowersetMap
pm_create(size_t bitmap_size /* In number of uint64_t elements */) {
    return (struct PowersetMap){.bitmap_size = bitmap_size,
                                .size = 0, .capacity = 512,
                                .keys = calloc(512, bitmap_size * sizeof(uint64_t)),
                                .states = calloc(512, sizeof(struct DFAState *))
    };
}

static void
pm_release(const struct PowersetMap *pm) {
    free(pm->keys);
    free(pm->states);
}

static size_t
pm_query(struct PowersetMap *pm, const uint64_t *key) {
    if (bitmap_iszero(key, pm->bitmap_size))
        return SIZE_MAX;
    for (size_t i = 0; i < pm->size; ++i) {
        if (bitmap_eq(key, pm->keys + i * pm->bitmap_size, pm->bitmap_size)) {
            return i;
        }
    }

    // Does not find preexisting
    // Make new
    if (pm->size >= pm->capacity) {
        size_t new_capacity = pm->capacity + 256;
        pm->keys = reallocarray(pm->keys, new_capacity, pm->bitmap_size * sizeof(uint64_t));
        pm->states = reallocarray(pm->states, new_capacity, sizeof(struct DFAState *));
        pm->capacity = new_capacity;
    }
    
    size_t value = pm->size++;
    pm->states[value] = NULL;
    memcpy(pm->keys + value * pm->bitmap_size, key, pm->bitmap_size * sizeof(uint64_t));
    // printf("QUERY MISS %lu for %lu%lu\n", value, *key, *(key+1));
    // printf("%lu%lu\n", test[0], test[1]);
    return value;
}

static size_t 
dfa_powerset_construction_dfs(const struct NFA *nfa, struct PowersetMap *pm, const uint64_t *bitmap) {
    size_t state_idx = pm_query(pm, bitmap);
    if (pm->states[state_idx] != NULL) {
        return state_idx;
    }
    pm->states[state_idx] = (struct DFAState *)1; // A place holder

    static_assert(sizeof(input) == sizeof(char), "Assume char");
    struct DFAState *state = dfa_state_init(NONTERMINAL, ASCII_TABLE_SIZE - 1);

    // Colate transitions
    size_t nondet_states[64],
           nondet_state_num = 0;
    for (size_t i = 0; i < nfa->num_state; ++i) {
        if (get_bit(bitmap, i)) {
            assert(nondet_state_num < 64); // If There are more than 64 nondeterminisitic states, then whole dfa would probably not fit in memory
            state->terminality |= nfa->states[i]->terminality;
            nondet_states[nondet_state_num++] = i;
        }
    }

    size_t transitions[64]; // Index where the nondet state transition are read

    for (size_t i = 0; i < nondet_state_num; ++i) {
        transitions[i] = 0;
        assert(nfa->states[nondet_states[i]]->num_transition == 0 || nfa->states[nondet_states[i]]->transitions[0].character != EPSILON);
    }

    uint64_t *private_bitmap = calloc(pm->bitmap_size, sizeof(uint64_t));
    struct Transition *write = state->unique_transitions;
    for (input ch = 1; isascii(ch); ++ch) {
        // Colate transitions
        for (size_t j = 0; j < nondet_state_num; ++j) {
            const struct State *state = nfa->states[nondet_states[j]];
            size_t k = transitions[j];
            for (; k < state->num_transition; ++k) {
                assert(state->transitions[k].character >= ch);
                if (state->transitions[k].character != ch)
                    break;
                set_bit(private_bitmap, state->transitions[k].successor, true);
            }
            transitions[j] = k;
        }

        // Write transitions
        if (!bitmap_iszero(private_bitmap, pm->bitmap_size))
            *(write++) = (struct Transition){.character = ch, .successor = dfa_powerset_construction_dfs(nfa, pm, private_bitmap)};
        
        // Reset bitmap
        memset(private_bitmap, 0, sizeof(uint64_t) * pm->bitmap_size);
    }
    free(private_bitmap);

    assert(write - state->unique_transitions <= state->num_transitions);
    state->num_transitions = write - state->unique_transitions;
    state = realloc(state, STRUCT_SIZE(struct DFAState, unique_transitions, state->num_transitions));

    // finalizing
    pm->states[state_idx] = state;
    return state_idx;
}

struct DFA *
dfa_powerset_construction(const struct NFA *nfa) {
    for (size_t i = 0; i < nfa->num_state; ++i)
        state_mut_sort_transition(nfa->states[i]);

    size_t bitmap_size = nfa->num_state / 64 + (nfa->num_state % 64 ? 1 : 0);
    struct PowersetMap pm = pm_create(bitmap_size);
    uint64_t *bitmap = calloc(bitmap_size, sizeof(uint64_t));
    for (size_t i = 0; i < nfa->num_initial; ++i)
        set_bit(bitmap, i, true);
    // printf("bitmap: %#lx\n", *bitmap);
    size_t start = dfa_powerset_construction_dfs(nfa, &pm, bitmap);
    free(bitmap);
    assert(start == 0);
    
    struct DFA *dfa = malloc(STRUCT_SIZE(struct DFA, states, pm.size));
    dfa->num_states = pm.size;
    memcpy(dfa->states, pm.states, sizeof(struct NFAState *) * pm.size);

    pm_release(&pm);
    return dfa;
}

void
dfa_free(struct DFA *nfa) {
    for (size_t i = 0; i < nfa->num_states; ++i)
        free(nfa->states[i]);
    free(nfa);
}

void 
dfa_fprint(FILE *stream, const struct DFA *automata) {
    fprintf(stream, "STARTING STATE: 0\n"
                    "TERMINAL STATE:\n");
    for (size_t i = 0; i < automata->num_states; ++i) {
        const struct DFAState *state = automata->states[i];
        if (state->terminality != NONTERMINAL) {
            fprintf(stream, "\tSTATE %#lx ", i);
            for (size_t i = 0; i < NUM_TOKEN_KIND; ++i)
                if (state->terminality & (1<<i)) 
                    fprintf(stream, "%s ", token_names[i]);
            fputc('\n', stream);
        }
    }
    fprintf(stream, "TRANSITIONS:\n");
    char cell[] = "                    ";
    fputs(cell, stream);
    for (input ch = 1; ch > 0; ++ch) {
        if (isprint(ch)) {
            cell[snprintf(cell, sizeof(cell) - 1, "'%c'", ch)] = ' ';
        } else {
            cell[snprintf(cell, sizeof(cell) - 1, "'\\%u'", ch)] = ' ';
        }
        fprintf(stream, "|%s", cell);
        memset(cell, ' ', sizeof(cell) - 1);
    }
    fputc('\n', stream);

    for (size_t i = 0; i < automata->num_states; ++i) {
        const struct DFAState *state = automata->states[i];
        cell[snprintf(cell, sizeof(cell) - 1, "%#lx", i)] = ' ';
        fputs(cell, stream);
        memset(cell, ' ', sizeof(cell) - 1);

        const struct Transition *trans = state->unique_transitions;
        for (input ch = 1; isascii(ch); ++ch) {
            assert(trans >= state->unique_transitions + state->num_transitions || ch <= trans->character);
            size_t successor = SIZE_MAX;
            if (trans < state->unique_transitions + state->num_transitions)
                successor = trans->character == ch ? (trans++)->successor : SIZE_MAX;
            
            if (successor != SIZE_MAX) {
                cell[snprintf(cell, sizeof(cell) - 1, "%#lx", successor)] = ' ';
            } else {
                cell[snprintf(cell, sizeof(cell) - 1, "NULL STATE")] = ' ';
            }

            fprintf(stream, "|%s", cell);
            memset(cell, ' ', sizeof(cell) - 1);
            
        }
        fputc('\n', stream);
    }
}

token 
dfa_recognize(const struct DFA *automata, input *str) {
    size_t state = 0;
    while (*str) {
        state = dfa_consume(automata, state, *str);
        ++str;
    }
    return dfa_state_token(automata, state);
}

size_t
dfa_consume(const struct DFA *automata, size_t state, input in) {
    if (state == SIZE_MAX) return state;
    assert(state < automata->num_states);
    return dfa_state_search(automata->states[state], in);
}

token
dfa_state_token(const struct DFA *automata, size_t state) {
    if (state == SIZE_MAX) return NONTERMINAL;
    assert(state < automata->num_states);
    return automata->states[state]->terminality;
}

static void
test1() {
    printf("TEST 1\n");
    struct NFA *test = nfa_verbatim("Hello, world!");
    nfa_fprint(stdout, test);
    nfa_free(test);
}

static void
test2() {
    printf("TEST 2\n");
    struct NFA *lhs = nfa_verbatim("Hello, world!"),
                               *rhs = nfa_verbatim("Goodbye, world!"),
                               *out = nfa_union(lhs, rhs);
    nfa_fprint(stdout, out);
    nfa_free(lhs);
    nfa_free(rhs);
    nfa_free(out);
}

static void 
test3() {
    printf("TEST 3\n");
    struct NFA *lhs = nfa_verbatim(""),
                               *rhs = nfa_verbatim("Goodbye, world!"),
                               *out = nfa_union(lhs, rhs);
    nfa_fprint(stdout, out);
    nfa_free(lhs);
    nfa_free(rhs);
    nfa_free(out);
}

static void
test4() {
    printf("TEST 4\n");
    struct NFA *prefix = nfa_verbatim("Hello, "),
                               *suffix = nfa_verbatim("World!"),
                               *output = nfa_concat(prefix, suffix);
    nfa_fprint(stdout, output);
    nfa_free(prefix);
    nfa_free(suffix);
    nfa_free(output);
}

static void
test5() {
    printf("TEST 5\n");
    struct NFA *prefix = nfa_verbatim(""),
                               *suffix = nfa_verbatim("World!"),
                               *output = nfa_concat(prefix, suffix);
    nfa_fprint(stdout, output);
    nfa_free(prefix);
    nfa_free(suffix);
    nfa_free(output);
}

static void
test6() {
    printf("TEST 6\n");
    struct NFA *input  = nfa_verbatim("Hello, World!"),
                               *output = nfa_kleene(input);
    nfa_fprint(stdout, output);
    nfa_free(input);
    nfa_free(output);
}

static void
test7() {
    printf("TEST 7\n");
    struct NFA *input  = nfa_verbatim("Hello, World!"),
               *output = nfa_kleene(input);
    nfa_mut_eliminate_epsilon(output);
    nfa_fprint(stdout, output);
    nfa_free(input);
    nfa_free(output);
}

static void
test9() {
    printf("TEST 9\n");
    struct NFA *alphanum = nfa_alphanum(),
               *lhs = nfa_verbatim("lhs"),
               *rhs = nfa_verbatim("rhs"),
               *a1  = nfa_kleene(alphanum),
               *a2  = nfa_concat(lhs, a1),
               *stf = nfa_verbatim("Hello, world!"),
               *a3  = nfa_concat(a2, rhs),
               *a4  = nfa_union(a3, stf);

    nfa_mut_minimize(a4);
    nfa_fprint(stdout, a4);
    // printf("%lu\n", a3->num_state);
    // nfa_fprint(stdout, a3);
    struct DFA *automata = dfa_powerset_construction(a4);
    FILE *output = fopen("testing", "w");
    dfa_fprint(output, automata);
    fclose(output);
    
#define TRY_AND_RECOGNIZE(str) printf(str " %lu\n", dfa_recognize(automata, str))
    TRY_AND_RECOGNIZE("Hello, world!");
    TRY_AND_RECOGNIZE("lhsrhs");
    TRY_AND_RECOGNIZE("lhsxrhs");
    TRY_AND_RECOGNIZE("lhs");

    TRY_AND_RECOGNIZE("egjywueyr");
#undef TRY_AND_RECOGINIZE

    dfa_free(automata);

    nfa_free(stf);
    nfa_free(a4);
    nfa_free(a3);
    nfa_free(a2);
    nfa_free(a1);
    nfa_free(rhs);
    nfa_free(lhs);
    nfa_free(alphanum);
}

