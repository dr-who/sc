/* rpn.h - Reverse Polish Notation stack engine header
 * C89 portable for DOS, Linux, VIC-20
 */

#ifndef RPN_H
#define RPN_H

#include "config.h"

#ifdef HAVE_RPN

#include "apfc.h"

/* Stack configuration */
#ifndef RPN_STACK_SIZE
#define RPN_STACK_SIZE 16
#endif

/* Stack storage - accessible externally */
extern apfc rpn_stack[RPN_STACK_SIZE];
extern int rpn_sp;
extern int rpn_mode;
extern apfc rpn_lastx;

/* Initialize RPN engine */
void rpn_init(void);

/* Stack operations */
void rpn_push(const apfc *val);
int rpn_pop(apfc *val);
int rpn_peek(apfc *val);
void rpn_drop(void);
void rpn_dup(void);
void rpn_swap(void);
void rpn_clear(void);
void rpn_roll_down(void);
void rpn_roll_up(void);

/* Display */
void rpn_show_stack(void);
void format_value_rpn(const apfc *val);

/* Process a single token (number, operator, or function) */
int rpn_process(const char *token);

/* Evaluate a full line of RPN input */
void rpn_eval_line(const char *line);

#endif /* HAVE_RPN */

#endif /* RPN_H */
