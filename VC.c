#include "VC.h"

#include "Lexer.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#define MAX_TOKEN_LENGTH 128
#define BUFFER_SIZE 512
#define MUNCH_LIMIT 512

// GLOBAL
struct DFA *vc_lexer_automata = NULL;

void
vc_lexer_automata_init() {
    // comment
    struct NFA *comment_start        = nfa_verbatim("/*"),
               *comment_end          = nfa_verbatim("*/"),
               *comment_typical_char = nfa_any(" !\"#$%&'()+,-.0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\n\t"),
               *comment_asterisk     = nfa_verbatim("*"),
               *comment_forwardslash = nfa_verbatim("/"),
               *comment_nonend       = nfa_concat(comment_asterisk, comment_typical_char),
               *comment_char         = nfa_union_many(comment_forwardslash, comment_typical_char, comment_nonend, NULL),
               *comment_char_star    = nfa_kleene(comment_char),
               *comment              = nfa_concat_many(comment_start, comment_char_star, comment_end, NULL);
    nfa_define_token(comment, TOKEN_IGNORE_BIT);

    // nonchar
    struct NFA *nonchar = nfa_any("\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10"
                                  "\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f\x7f "),
               *nonchar_star = nfa_kleene(nonchar),
               *nonchar_ = nfa_concat(nonchar, nonchar_star);
    nfa_define_token(nonchar_, TOKEN_IGNORE_BIT);

    // atoms
    struct NFA *letter    = nfa_any(LETTER),
               *digit     = nfa_any(DIGIT),
               *alphanum  = nfa_any(DIGIT LETTER),
               *seperator = nfa_any("{}()[];,");
    nfa_define_token(seperator, TOKEN_SEPERATOR_BIT);
    
    // kleene
    struct NFA *alphanum_star = nfa_kleene(alphanum),
               *identifier    = nfa_concat(letter, alphanum_star);
    nfa_define_token(identifier, TOKEN_IDENTIFIER_BIT);
//    nfa_fprint(stdout, identifier);

    // keywords
    struct NFA *keyword_boolean  = nfa_verbatim("boolean"),
               *keyword_break    = nfa_verbatim("break"),
               *keyword_continue = nfa_verbatim("continue"),
               *keyword_else     = nfa_verbatim("else"),
               *keyword_for      = nfa_verbatim("for"),
               *keyword_float    = nfa_verbatim("float"),
               *keyword_if       = nfa_verbatim("if"),
               *keyword_int      = nfa_verbatim("int"),
               *keyword_return   = nfa_verbatim("return"),
               *keyword_void     = nfa_verbatim("void"),
               *keyword_while    = nfa_verbatim("while"),
               *keyword          = nfa_union_many(keyword_boolean, keyword_break, keyword_continue, keyword_else, keyword_for, keyword_float, keyword_if, keyword_int, keyword_return, keyword_void, keyword_while, NULL);
    nfa_define_token(keyword, TOKEN_KEYWORD_BIT);

    // Operator
    struct NFA *operator_singlechar = nfa_any("+-*/=!><"),
               *operator_ge         = nfa_verbatim(">="),
               *operator_le         = nfa_verbatim("<="),
               *operator_ne         = nfa_verbatim("!="),
               *operator_eq         = nfa_verbatim("=="),
               *operator_or         = nfa_verbatim("||"),
               *operator_and        = nfa_verbatim("&&"),
               *operator = nfa_union_many(operator_singlechar, operator_ge, operator_le, operator_ne, operator_eq, operator_or, operator_and, NULL);
    nfa_define_token(operator, TOKEN_OPERATOR_BIT);

    // Literals
    struct NFA *digit_star = nfa_kleene(digit),
               *point      = nfa_verbatim("."),
               *point_     = nfa_optional(point),
               *literal_int = nfa_concat(digit, digit_star),
               *exponent_e = nfa_any("eE"),
               *exponent_s = nfa_any("+-"),
               *exponent_s_ = nfa_optional(exponent_s),
               *exponent    = nfa_concat_many(exponent_e, exponent_s_, literal_int, NULL), 
               *exponent_   = nfa_optional(exponent),
               *fraction    = nfa_concat(point, literal_int),
               *literal_float1 = nfa_concat_many(digit_star, fraction, exponent_, NULL),
               *literal_float2 = nfa_concat(literal_int, point),
               *literal_float3 = nfa_concat_many(literal_int, point_, exponent, NULL),
               *literal_float  = nfa_union_many(literal_float1, literal_float2, literal_float3, NULL),
               *literal_true  = nfa_verbatim("true"),
               *literal_false = nfa_verbatim("false"),
               *literal_boolean = nfa_union(literal_true, literal_false);
    struct NFA *backslash  = nfa_verbatim("\\"), 
               *escapeable = nfa_any("\\bfnrt'\""),
               *nonescape = nfa_any(DIGIT LETTER " !#$%&'()*+,-./:;<>=?@[]^_`~"),
               *escaped = nfa_concat(backslash, escapeable),
               *string_interal = nfa_union(escaped, nonescape),
               *string_interal_star = nfa_kleene(string_interal),
               *quote = nfa_verbatim("\""),
               *literal_string = nfa_concat_many(quote, string_interal_star, quote, NULL);
    struct NFA *literal = nfa_union_many(literal_int, literal_float, literal_boolean, literal_string, NULL);
    nfa_define_token(literal, TOKEN_LITERAL_BIT);

    struct NFA *nfa = nfa_union_many(comment, identifier, keyword, operator, seperator, literal, nonchar_, NULL);
    nfa_mut_minimize(nfa);
 //   nfa_fprint(stdout, nfa);

    vc_lexer_automata = dfa_powerset_construction(nfa);
    nfa_free(nfa);

    // Literals
    nfa_free(literal);

    nfa_free(literal_string);
    nfa_free(quote);
    nfa_free(string_interal_star);
    nfa_free(string_interal);
    nfa_free(nonescape);
    nfa_free(escaped);
    nfa_free(escapeable);
    nfa_free(backslash);

    nfa_free(literal_boolean);
    nfa_free(literal_true);
    nfa_free(literal_false);

    nfa_free(literal_float);

    nfa_free(literal_float1);
    nfa_free(literal_float2);
    nfa_free(literal_float3);

    nfa_free(exponent_);
    nfa_free(point_);
    nfa_free(fraction);

    nfa_free(exponent_e);
    nfa_free(exponent_s);
    nfa_free(exponent_s_);
    nfa_free(exponent);

    nfa_free(literal_int);
    nfa_free(digit_star);
    nfa_free(point);

    // Operator
    nfa_free(operator);
    nfa_free(operator_singlechar);
    nfa_free(operator_ge        );
    nfa_free(operator_le        );
    nfa_free(operator_ne        );
    nfa_free(operator_eq        );
    nfa_free(operator_or        );
    nfa_free(operator_and       );





    // keywords
    nfa_free(keyword);

    nfa_free(keyword_boolean );
    nfa_free(keyword_break   );
    nfa_free(keyword_continue);
    nfa_free(keyword_else    );
    nfa_free(keyword_for     );
    nfa_free(keyword_float   );
    nfa_free(keyword_if      );
    nfa_free(keyword_int     );
    nfa_free(keyword_return  );
    nfa_free(keyword_void    );
    nfa_free(keyword_while   );

    
    // kleene
    nfa_free(identifier);
    nfa_free(alphanum_star);

    // atoms
    nfa_free(seperator);
    nfa_free(alphanum);
    nfa_free(letter);
    nfa_free(digit);

    // nonchar
    nfa_free(nonchar_);
    nfa_free(nonchar_star);
    nfa_free(nonchar);

    // comment
    nfa_free(comment_start        );
    nfa_free(comment_end          );
    nfa_free(comment_typical_char );
    nfa_free(comment_asterisk     );
    nfa_free(comment_forwardslash );
    nfa_free(comment_nonend       );
    nfa_free(comment_char         );
    nfa_free(comment_char_star    );
    nfa_free(comment              );
}


void scanfile(FILE *in, FILE *out) {
    input munch[MUNCH_LIMIT];
    size_t state = 0,
           token_length = 0;  // Length of longest known token
    token  tok = NONTERMINAL; // Kind of the longest known token
    long   start = ftell(in);
    assert(dfa_state_token(vc_lexer_automata, 0) == NONTERMINAL);
    
    size_t nread = fread(munch, sizeof(input), MUNCH_LIMIT, in);
    if (nread == 0)
        return;
    
    size_t consumed = 0;
    bool found_token = false;
    for (const input *ch = munch; ch < munch + nread; ++ch) {
        ++consumed;
        size_t old_state = state;
        state = dfa_consume(vc_lexer_automata, state, *ch);
        if (state == SIZE_MAX) {
            found_token = true;
            if (tok == NONTERMINAL) { // TODO: Add spaces to ends of tokens so remove this hack TODO: Remove this
                fprintf(stderr, "UNIDENTIFIABLE TOKEN: \"");
                fwrite(munch, sizeof(input), consumed, stderr);
                fputs("\"\n", stderr);
                fprintf(stderr, "STATE: %#0lx CHAR: %#0x\n", old_state, *ch);
                exit(-1);
            } else if (tok != TOKEN_IGNORE_BIT) {
                assert(TOKEN_IGNORE_BIT & ~tok); // Must not have this bit set
                fwrite(munch, sizeof(input), token_length, out);
                fputc('\t', out);
                for (size_t i = 0; i < NUM_TOKEN_KIND; ++i) 
                    if (tok & (1<<i))
                        fprintf(out, "%s ", token_names[i]);
                fputc('\n', out);
            }
            break;
        }

        token new_tok = dfa_state_token(vc_lexer_automata, state);
        if (new_tok != NONTERMINAL) {
            token_length = consumed;
            tok = new_tok;
        }
    }

    // handle EOF
    if (!found_token && feof(in)) {
        if (tok == NONTERMINAL) {
            fprintf(stderr, "UNIDENTIFIABLE TOKEN: ");
            fwrite(munch, sizeof(input), token_length, stderr);
            fputc('\n', stderr);
            exit(-1);
        } else if (tok != TOKEN_IGNORE_BIT) {
            assert(TOKEN_IGNORE_BIT & ~tok);
            fwrite(munch, sizeof(input), token_length, out);
            fputc('\t', out);
            for (size_t i = 0; i < NUM_TOKEN_KIND; ++i) 
                if (tok & (1<<i))
                    fprintf(out, "%s ", token_names[i]);
            fputc('\n', out);
        }
        return;
    }

    assert(found_token);
    fseek(in, start + token_length, SEEK_SET);
    scanfile(in, out);
}

int main(int argc, char *argv[argc]) {
    if (argc < 3) {
        printf("USAGE: %s <*.vc> <*.vctok>\n", argv[0]);
        return 0;
    }
    vc_lexer_automata_init();
    
    FILE *source = fopen(argv[1], "r");
    FILE *dest   = fopen(argv[2], "w");
    scanfile(source, dest);
    assert(feof(source));
    fclose(source);
    fclose(dest);
    vc_lexer_automata_free();
    return 0;
}
