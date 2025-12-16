/*
 * func_dispatch.c - Function documentation lookup and dispatch system
 * 
 * Provides unified access to all function documentation tables.
 * Used by help, demo, test commands, and documentation generation.
 *
 * C89 compliant
 */
#include "func_defs.h"
#include <stdio.h>
#include <string.h>

/* Import all function tables */
extern const FuncDoc func_trig_docs[];
extern const FuncDoc func_exp_docs[];
extern const FuncDoc func_log_docs[];
extern const FuncDoc func_complex_docs[];
extern const FuncDoc func_round_docs[];
extern const FuncDoc func_number_docs[];
extern const FuncDoc func_special_docs[];
extern const FuncDoc func_stats_docs[];
extern const FuncDoc func_prob_docs[];
extern const FuncDoc func_linalg_docs[];
extern const FuncDoc func_matrix_docs[];
extern const FuncDoc func_matquery_docs[];
extern const FuncDoc func_signal_docs[];
extern const FuncDoc func_poly_docs[];
extern const FuncDoc func_angle_docs[];
extern const FuncDoc func_set_docs[];
extern const FuncDoc func_logic_docs[];
extern const FuncDoc func_bitwise_docs[];
extern const FuncDoc func_coord_docs[];
extern const FuncDoc func_util_docs[];
extern const FuncDoc func_const_docs[];
extern const FuncDoc func_datasci_docs[];
extern const FuncDoc func_text_docs[];

/* All doc tables for searching */
static const FuncDoc *all_tables[] = {
    NULL, /* func_trig_docs - set in init */
    NULL, /* func_exp_docs */
    NULL, /* func_log_docs */
    NULL, /* func_complex_docs */
    NULL, /* func_round_docs */
    NULL, /* func_number_docs */
    NULL, /* func_special_docs */
    NULL, /* func_stats_docs */
    NULL, /* func_prob_docs */
    NULL, /* func_linalg_docs */
    NULL, /* func_matrix_docs */
    NULL, /* func_matquery_docs */
    NULL, /* func_signal_docs */
    NULL, /* func_poly_docs */
    NULL, /* func_angle_docs */
    NULL, /* func_set_docs */
    NULL, /* func_logic_docs */
    NULL, /* func_bitwise_docs */
    NULL, /* func_coord_docs */
    NULL, /* func_util_docs */
    NULL, /* func_const_docs */
    NULL, /* func_datasci_docs */
    NULL, /* func_text_docs */
    NULL  /* Sentinel */
};

/* Initialization flag */
static int tables_initialized = 0;

/* Initialize all table pointers */
static void init_tables(void)
{
    if (tables_initialized) return;
    
    all_tables[0] = func_trig_docs;
    all_tables[1] = func_exp_docs;
    all_tables[2] = func_log_docs;
    all_tables[3] = func_complex_docs;
    all_tables[4] = func_round_docs;
    all_tables[5] = func_number_docs;
    all_tables[6] = func_special_docs;
    all_tables[7] = func_stats_docs;
    all_tables[8] = func_prob_docs;
    all_tables[9] = func_linalg_docs;
    all_tables[10] = func_matrix_docs;
    all_tables[11] = func_matquery_docs;
    all_tables[12] = func_signal_docs;
    all_tables[13] = func_poly_docs;
    all_tables[14] = func_angle_docs;
    all_tables[15] = func_set_docs;
    all_tables[16] = func_logic_docs;
    all_tables[17] = func_bitwise_docs;
    all_tables[18] = func_coord_docs;
    all_tables[19] = func_util_docs;
    all_tables[20] = func_const_docs;
    all_tables[21] = func_datasci_docs;
    all_tables[22] = func_text_docs;
    all_tables[23] = NULL;
    
    tables_initialized = 1;
}

/* ========== Lookup Functions ========== */

const FuncDoc *func_find_doc(const char *name)
{
    int t, i;
    
    init_tables();
    
    for (t = 0; all_tables[t] != NULL; t++) {
        const FuncDoc *table = all_tables[t];
        for (i = 0; table[i].name != NULL; i++) {
            if (strcmp(table[i].name, name) == 0) {
                return &table[i];
            }
        }
    }
    return NULL;
}

int func_doc_exists(const char *name)
{
    return func_find_doc(name) != NULL;
}

int func_doc_count(void)
{
    int t, i;
    int count = 0;
    
    init_tables();
    
    for (t = 0; all_tables[t] != NULL; t++) {
        const FuncDoc *table = all_tables[t];
        for (i = 0; table[i].name != NULL; i++) {
            count++;
        }
    }
    return count;
}

/* ========== Help System ========== */

void func_show_help(const char *name)
{
    const FuncDoc *doc = func_find_doc(name);
    int i;
    
    if (doc == NULL) {
        printf("No help available for '%s'\n", name);
        printf("Type 'functions' to list all available functions.\n");
        return;
    }
    
    printf("\n%s - %s\n", doc->name, doc->category);
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("Syntax:      %s\n", doc->syntax);
    printf("Description: %s\n", doc->description);
    
    printf("\nExamples:\n");
    for (i = 0; i < MAX_EXAMPLES && doc->examples[i] != NULL; i++) {
        printf("  %s\n", doc->examples[i]);
    }
    
    if (doc->see_also != NULL && strlen(doc->see_also) > 0) {
        printf("\nSee also: %s\n", doc->see_also);
    }
    printf("\n");
}

void func_show_demo(const char *name)
{
    const FuncDoc *doc = func_find_doc(name);
    int i;
    extern int eval_expr_line(const char *line, int quiet);
    
    if (doc == NULL) {
        printf("No demo available for '%s'\n", name);
        return;
    }
    
    printf("\n=== Demo: %s ===\n", doc->name);
    printf("Description: %s\n\n", doc->description);
    
    for (i = 0; i < MAX_EXAMPLES && doc->examples[i] != NULL; i++) {
        const char *ex = doc->examples[i];
        const char *eq = strchr(ex, '=');
        char expr[256];
        
        /* Extract just the expression (before =) */
        if (eq != NULL) {
            int len = (int)(eq - ex);
            if (len > 0 && len < 255) {
                strncpy(expr, ex, len);
                expr[len] = '\0';
                /* Trim trailing space */
                while (len > 0 && expr[len-1] == ' ') {
                    expr[--len] = '\0';
                }
            } else {
                strncpy(expr, ex, 255);
                expr[255] = '\0';
            }
        } else {
            strncpy(expr, ex, 255);
            expr[255] = '\0';
        }
        
        printf(">>> %s\n", expr);
        eval_expr_line(expr, 0);
        printf("\n");
    }
}

/* ========== List Functions ========== */

void func_list_all(void)
{
    const char *printed_cats[30];
    int num_printed = 0;
    int t, i, j, c;
    
    init_tables();
    
    printf("\nBuilt-in Functions (%d total):\n", func_doc_count());
    printf("=================================\n\n");
    
    /* Print by category */
    for (t = 0; all_tables[t] != NULL; t++) {
        const FuncDoc *table = all_tables[t];
        
        for (i = 0; table[i].name != NULL; i++) {
            const char *cat = table[i].category;
            int already = 0;
            
            /* Check if we've printed this category */
            for (c = 0; c < num_printed; c++) {
                if (strcmp(printed_cats[c], cat) == 0) {
                    already = 1;
                    break;
                }
            }
            
            if (!already && num_printed < 30) {
                int count = 0;
                
                printf("%s:\n  ", cat);
                printed_cats[num_printed++] = cat;
                
                /* Print all functions in this category */
                for (j = 0; all_tables[j] != NULL; j++) {
                    const FuncDoc *t2 = all_tables[j];
                    int k;
                    for (k = 0; t2[k].name != NULL; k++) {
                        if (strcmp(t2[k].category, cat) == 0) {
                            printf("%s ", t2[k].name);
                            count++;
                            if (count % 10 == 0) {
                                printf("\n  ");
                            }
                        }
                    }
                }
                printf("\n\n");
            }
        }
    }
    
    printf("Use 'help <func>' for details, 'demo <func>' for examples.\n\n");
}

/* ========== Test Generation ========== */

int func_run_tests(void)
{
    int t, i, j;
    int passed = 0, failed = 0;
    extern int eval_expr_line(const char *line, int quiet);
    
    init_tables();
    
    printf("\nRunning Function Tests from Documentation:\n");
    printf("==========================================\n\n");
    
    for (t = 0; all_tables[t] != NULL; t++) {
        const FuncDoc *table = all_tables[t];
        
        for (i = 0; table[i].name != NULL; i++) {
            const FuncDoc *doc = &table[i];
            
            for (j = 0; j < MAX_EXAMPLES && doc->examples[j] != NULL; j++) {
                const char *ex = doc->examples[j];
                const char *eq;
                char expr[256];
                int result;
                
                /* Skip examples that are just showing syntax (no =) */
                eq = strchr(ex, '=');
                if (eq == NULL) continue;
                
                /* Extract expression part */
                {
                    int len = (int)(eq - ex);
                    if (len > 0 && len < 255) {
                        strncpy(expr, ex, len);
                        expr[len] = '\0';
                        while (len > 0 && expr[len-1] == ' ') {
                            expr[--len] = '\0';
                        }
                    } else {
                        continue;
                    }
                }
                
                /* Skip matrix/multi-line examples */
                if (strchr(expr, ';') != NULL && strchr(expr, '[') == NULL) continue;
                if (strstr(expr, "= eye") != NULL) continue;
                if (strstr(expr, "= A") != NULL) continue;
                
                /* Run the expression */
                result = eval_expr_line(expr, 1);
                if (result >= 0) {
                    passed++;
                } else {
                    printf("  FAIL: %s\n", expr);
                    failed++;
                }
            }
        }
    }
    
    printf("\nDocumentation tests: %d passed, %d failed\n\n", passed, failed);
    return failed;
}

/* ========== Documentation Generation ========== */

void func_gen_test_script(void)
{
    int t, i, j;
    
    init_tables();
    
    printf("#!/bin/bash\n");
    printf("# Auto-generated from function documentation\n");
    printf("# Run: bash tests/test_funcs_gen.sh\n\n");
    printf("SC=\"${SC:-./bin/sc}\"\n");
    printf("PASS=0\nFAIL=0\n\n");
    printf("test_expr() {\n");
    printf("    result=$($SC -e \"$1\" 2>&1)\n");
    printf("    if [ $? -eq 0 ]; then\n");
    printf("        echo \"  PASS: $1\"\n");
    printf("        PASS=$((PASS + 1))\n");
    printf("    else\n");
    printf("        echo \"  FAIL: $1\"\n");
    printf("        FAIL=$((FAIL + 1))\n");
    printf("    fi\n");
    printf("}\n\n");
    
    for (t = 0; all_tables[t] != NULL; t++) {
        const FuncDoc *table = all_tables[t];
        
        for (i = 0; table[i].name != NULL; i++) {
            const FuncDoc *doc = &table[i];
            
            for (j = 0; j < MAX_EXAMPLES && doc->examples[j] != NULL; j++) {
                const char *ex = doc->examples[j];
                const char *eq = strchr(ex, '=');
                
                if (eq != NULL) {
                    char expr[256];
                    int len = (int)(eq - ex);
                    if (len > 0 && len < 255) {
                        strncpy(expr, ex, len);
                        expr[len] = '\0';
                        while (len > 0 && expr[len-1] == ' ') {
                            expr[--len] = '\0';
                        }
                        /* Escape special chars for shell */
                        printf("test_expr \"%s\"\n", expr);
                    }
                }
            }
        }
    }
    
    printf("\necho \"\"\n");
    printf("echo \"Results: $PASS passed, $FAIL failed\"\n");
    printf("exit $FAIL\n");
}

void func_gen_manpage(void)
{
    const char *printed_cats[30];
    int num_printed = 0;
    int t, i, j, c;
    
    init_tables();
    
    printf(".\\\" Auto-generated function reference\n");
    printf(".SH FUNCTIONS\n");
    printf("This section documents all available functions organized by category.\n");
    
    for (t = 0; all_tables[t] != NULL; t++) {
        const FuncDoc *table = all_tables[t];
        
        for (i = 0; table[i].name != NULL; i++) {
            const char *cat = table[i].category;
            int already = 0;
            
            for (c = 0; c < num_printed; c++) {
                if (strcmp(printed_cats[c], cat) == 0) {
                    already = 1;
                    break;
                }
            }
            
            if (!already && num_printed < 30) {
                printf(".SS %s\n", cat);
                printed_cats[num_printed++] = cat;
                
                for (j = 0; all_tables[j] != NULL; j++) {
                    const FuncDoc *t2 = all_tables[j];
                    int k;
                    for (k = 0; t2[k].name != NULL; k++) {
                        if (strcmp(t2[k].category, cat) == 0) {
                            const FuncDoc *doc = &t2[k];
                            printf(".TP\n");
                            printf(".B %s\n", doc->syntax);
                            printf("%s\n", doc->description);
                        }
                    }
                }
            }
        }
    }
}
