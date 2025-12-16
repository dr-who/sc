/*
 * func_registry.c - Function registration and documentation system
 * C89 compliant
 */
#include "func_registry.h"
#include <stdio.h>
#include <string.h>

/* External category declarations */
extern const func_category_t func_linalg_category;
extern const func_category_t func_datasci_category;
/* Add more categories here as they're created:
extern const func_category_t func_trig_category;
extern const func_category_t func_stats_category;
extern const func_category_t func_matrix_category;
extern const func_category_t func_special_category;
*/

/* Global category list */
static const func_category_t *all_categories[] = {
    &func_linalg_category,
    &func_datasci_category,
    /* Add more categories here */
    NULL
};

/* Count of categories */
static int category_count = 0;

/* Initialize the registry */
void func_init_all(void)
{
    int i;
    category_count = 0;
    for (i = 0; all_categories[i] != NULL; i++) {
        category_count++;
    }
}

/* Find a function by name across all categories */
const func_entry_t *func_find(const char *name)
{
    int i, j;
    for (i = 0; all_categories[i] != NULL; i++) {
        const func_category_t *cat = all_categories[i];
        for (j = 0; j < cat->count; j++) {
            if (strcmp(cat->functions[j].name, name) == 0) {
                return &cat->functions[j];
            }
        }
    }
    return NULL;
}

/* Find function in a specific category */
const func_entry_t *func_find_in_category(const char *category, const char *name)
{
    int i, j;
    for (i = 0; all_categories[i] != NULL; i++) {
        const func_category_t *cat = all_categories[i];
        if (strcmp(cat->name, category) == 0) {
            for (j = 0; j < cat->count; j++) {
                if (strcmp(cat->functions[j].name, name) == 0) {
                    return &cat->functions[j];
                }
            }
            return NULL;
        }
    }
    return NULL;
}

/* List all functions (for 'functions' command) */
void func_list_all(void)
{
    int i, j;
    int total = 0;
    
    printf("\nRegistered Functions by Category:\n");
    printf("==================================\n\n");
    
    for (i = 0; all_categories[i] != NULL; i++) {
        const func_category_t *cat = all_categories[i];
        printf("%s (%d functions):\n", cat->name, cat->count);
        printf("  ");
        for (j = 0; j < cat->count; j++) {
            printf("%s", cat->functions[j].name);
            if (j < cat->count - 1) printf(", ");
            if ((j + 1) % 8 == 0 && j < cat->count - 1) printf("\n  ");
        }
        printf("\n\n");
        total += cat->count;
    }
    
    printf("Total: %d functions in %d categories\n\n", total, category_count);
    printf("Use 'help <function>' for details, 'demo <function>' for examples.\n\n");
}

/* List functions in a specific category */
void func_list_category(const char *category)
{
    int i, j;
    for (i = 0; all_categories[i] != NULL; i++) {
        const func_category_t *cat = all_categories[i];
        if (strcmp(cat->name, category) == 0) {
            printf("\n%s Functions:\n", cat->name);
            for (j = 0; j < cat->count; j++) {
                printf("  %-15s %s\n", 
                       cat->functions[j].syntax,
                       cat->functions[j].description);
            }
            printf("\n");
            return;
        }
    }
    printf("Unknown category: %s\n", category);
}

/* Count total functions */
int func_count_all(void)
{
    int i;
    int total = 0;
    for (i = 0; all_categories[i] != NULL; i++) {
        total += all_categories[i]->count;
    }
    return total;
}

/* Show help for a function */
void func_show_help(const char *name)
{
    const func_entry_t *f = func_find(name);
    int i;
    
    if (f == NULL) {
        printf("Unknown function: %s\n", name);
        printf("Use 'functions' to list all available functions.\n");
        return;
    }
    
    printf("\n%s - %s\n", f->name, f->category);
    printf("----------------------------------------\n");
    printf("Syntax:      %s\n", f->syntax);
    printf("Description: %s\n", f->description);
    
    if (f->examples[0] != NULL) {
        printf("\nExamples:\n");
        for (i = 0; i < 4 && f->examples[i] != NULL; i++) {
            printf("  %s\n", f->examples[i]);
        }
    }
    printf("\n");
}

/* Show demo for a function (runs examples) */
void func_show_demo(const char *name)
{
    const func_entry_t *f = func_find(name);
    int i;
    extern int eval_expr_line(const char *line, int quiet);
    
    if (f == NULL) {
        printf("Unknown function: %s\n", name);
        return;
    }
    
    printf("\n=== Demo: %s ===\n", f->name);
    printf("Description: %s\n\n", f->description);
    
    for (i = 0; i < 4 && f->examples[i] != NULL; i++) {
        printf(">>> %s\n", f->examples[i]);
        eval_expr_line(f->examples[i], 0);
        printf("\n");
    }
}

/* Run tests on all registered functions */
int func_run_tests(void)
{
    int i, j, k;
    int passed = 0, failed = 0;
    extern int eval_expr_line(const char *line, int quiet);
    
    printf("\nTesting Registered Functions:\n");
    printf("=============================\n\n");
    
    for (i = 0; all_categories[i] != NULL; i++) {
        const func_category_t *cat = all_categories[i];
        printf("Testing %s...\n", cat->name);
        
        for (j = 0; j < cat->count; j++) {
            const func_entry_t *f = &cat->functions[j];
            int func_ok = 1;
            
            /* Run each example */
            for (k = 0; k < 4 && f->examples[k] != NULL; k++) {
                int result = eval_expr_line(f->examples[k], 1);
                if (result != 0) {
                    printf("  FAIL: %s: %s\n", f->name, f->examples[k]);
                    func_ok = 0;
                }
            }
            
            if (func_ok) {
                passed++;
            } else {
                failed++;
            }
        }
    }
    
    printf("\nResults: %d passed, %d failed\n\n", passed, failed);
    return failed;
}

/* Benchmark all registered functions */
void func_run_bench(void)
{
    int i;
    int total = func_count_all();
    
    printf("\nFunction Benchmark\n");
    printf("==================\n\n");
    printf("Total registered functions: %d\n", total);
    printf("Categories: %d\n\n", category_count);
    
    for (i = 0; all_categories[i] != NULL; i++) {
        printf("  %s: %d functions\n", 
               all_categories[i]->name, 
               all_categories[i]->count);
    }
    printf("\n");
    
    /* TODO: Add timing benchmarks */
}
