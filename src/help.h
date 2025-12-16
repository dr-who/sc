/*
 * help.h - Function help and demo system
 */
#ifndef HELP_H
#define HELP_H

#include "func_defs.h"

/* Show help for a specific function */
void show_function_help(const char *name);

/* Show demo for a specific function */
void show_function_demo(const char *name);

/* Run benchmark on all functions */
void run_function_bench(void);

/* Check if function exists */
int function_exists(const char *name);

/* List all functions by category */
void list_all_functions(void);

/* Count total documented functions */
int count_help_entries(void);

#endif /* HELP_H */
