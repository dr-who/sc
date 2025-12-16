/* repl.c - REPL loop, line editing, and history
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"

/* ========== DOS Screen Setup ========== */

#ifdef HAVE_CONIO
/* DOS: Set up blue background with yellow text, clear screen */
void dos_init_screen(void)
{
    union REGS regs;
    unsigned int far *vidmem;
    int i;
    
    /* Clear screen using INT 10h, AH=06h (scroll up) */
    regs.h.ah = 0x06;       /* Scroll up */
    regs.h.al = 0x00;       /* Clear entire window */
    regs.h.bh = 0x1E;       /* Attribute: blue bg (1) + yellow fg (E) */
    regs.h.ch = 0;          /* Upper left row */
    regs.h.cl = 0;          /* Upper left col */
    regs.h.dh = 24;         /* Lower right row */
    regs.h.dl = 79;         /* Lower right col */
    int86(0x10, &regs, &regs);
    
    /* Move cursor to top-left */
    regs.h.ah = 0x02;       /* Set cursor position */
    regs.h.bh = 0;          /* Page 0 */
    regs.h.dh = 0;          /* Row 0 */
    regs.h.dl = 0;          /* Col 0 */
    int86(0x10, &regs, &regs);
    
    /* Set text attribute for all subsequent output */
    /* Write spaces with attribute to fill screen */
    vidmem = (unsigned int far *)MK_FP(0xB800, 0);
    for (i = 0; i < 2000; i++) {  /* 80x25 = 2000 cells */
        vidmem[i] = 0x1E20;  /* Yellow on blue, space char */
    }
    
    /* Reset cursor to top */
    regs.h.ah = 0x02;
    regs.h.bh = 0;
    regs.h.dh = 0;
    regs.h.dl = 0;
    int86(0x10, &regs, &regs);
}

void print_prompt(void)
{
    printf(">>> ");
}
#define HAVE_COLOR_PROMPT 0

#elif defined(HAVE_TERMIOS)
/* Unix/Linux: Use ANSI escape codes */
void dos_init_screen(void)
{
    /* No-op on Linux */
}

void print_prompt(void)
{
    printf("\033[33m>>>\033[0m ");
    fflush(stdout);
}
#define HAVE_COLOR_PROMPT 1

#else
/* Fallback */
void dos_init_screen(void)
{
}

void print_prompt(void)
{
    printf(">>> ");
}
#define HAVE_COLOR_PROMPT 0
#endif

/* ========== Line Editing ========== */

static char history[MAX_HISTORY][MAX_INPUT];
static int history_count = 0;
static int history_pos = 0;

#ifdef HAVE_TERMIOS
static struct termios orig_termios;
static int raw_mode = 0;

static void disable_raw_mode(void)
{
    if (raw_mode) {
        tcsetattr(STDIN_FILENO, TCSAFLUSH, &orig_termios);
        raw_mode = 0;
    }
}

static void enable_raw_mode(void)
{
    struct termios raw;
    if (!raw_mode) {
        tcgetattr(STDIN_FILENO, &orig_termios);
        raw = orig_termios;
        raw.c_lflag &= ~(ICANON | ECHO);
        raw.c_cc[VMIN] = 1;
        raw.c_cc[VTIME] = 0;
        tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
        raw_mode = 1;
    }
}

int read_line(char *buf, int max_len)
{
    int pos = 0;
    int len = 0;
    int c;
    
    /* For pipe/redirect input, use simple fgets */
    if (!isatty(STDIN_FILENO)) {
        if (fgets(buf, max_len, stdin) == NULL) {
            return -1;
        }
        len = strlen(buf);
        if (len > 0 && buf[len-1] == '\n') {
            buf[--len] = '\0';
        }
        return len;
    }
    
    enable_raw_mode();
    history_pos = history_count;
    
    while (1) {
        c = getchar();
        
        if (c == EOF || c == 4) {  /* Ctrl-D */
            disable_raw_mode();
            return -1;
        }
        
        if (c == '\n' || c == '\r') {
            printf("\n");
            buf[len] = '\0';
            
            /* Add to history if non-empty */
            if (len > 0 && (history_count == 0 || strcmp(buf, history[history_count - 1]) != 0)) {
                if (history_count < MAX_HISTORY) {
                    strcpy(history[history_count++], buf);
                } else {
                    int i;
                    for (i = 0; i < MAX_HISTORY - 1; i++) {
                        strcpy(history[i], history[i + 1]);
                    }
                    strcpy(history[MAX_HISTORY - 1], buf);
                }
            }
            
            disable_raw_mode();
            return len;
        }
        
        if (c == 127 || c == 8) {  /* Backspace */
            if (pos > 0) {
                int i;
                for (i = pos - 1; i < len - 1; i++) {
                    buf[i] = buf[i + 1];
                }
                len--;
                pos--;
                
                /* Redraw */
                printf("\r\033[33m>>>\033[0m ");
                buf[len] = '\0';
                printf("%s ", buf);
                printf("\r\033[33m>>>\033[0m ");
                for (i = 0; i < pos; i++) putchar(buf[i]);
                for (i = pos; i < len; i++) putchar(buf[i]);
                for (i = len; i < len + 1; i++) putchar(' ');
                printf("\r\033[33m>>>\033[0m ");
                for (i = 0; i < pos; i++) putchar(buf[i]);
            }
            continue;
        }
        
        if (c == 27) {  /* Escape sequence */
            int c2 = getchar();
            if (c2 == '[') {
                int c3 = getchar();
                if (c3 == 'A') {  /* Up arrow */
                    if (history_pos > 0) {
                        history_pos--;
                        strcpy(buf, history[history_pos]);
                        len = strlen(buf);
                        pos = len;
                        printf("\r\033[33m>>>\033[0m \033[K%s", buf);
                    }
                } else if (c3 == 'B') {  /* Down arrow */
                    if (history_pos < history_count - 1) {
                        history_pos++;
                        strcpy(buf, history[history_pos]);
                        len = strlen(buf);
                        pos = len;
                        printf("\r\033[33m>>>\033[0m \033[K%s", buf);
                    } else if (history_pos == history_count - 1) {
                        history_pos = history_count;
                        buf[0] = '\0';
                        len = 0;
                        pos = 0;
                        printf("\r\033[33m>>>\033[0m \033[K");
                    }
                } else if (c3 == 'C') {  /* Right arrow */
                    if (pos < len) {
                        putchar(buf[pos]);
                        pos++;
                    }
                } else if (c3 == 'D') {  /* Left arrow */
                    if (pos > 0) {
                        pos--;
                        printf("\033[D");
                    }
                }
            }
            continue;
        }
        
        if (c >= 32 && c < 127 && len < max_len - 1) {
            int i;
            for (i = len; i > pos; i--) {
                buf[i] = buf[i - 1];
            }
            buf[pos] = c;
            len++;
            pos++;
            
            buf[len] = '\0';
            printf("\r\033[33m>>>\033[0m %s", buf);
            printf("\r\033[33m>>>\033[0m ");
            for (i = 0; i < pos; i++) putchar(buf[i]);
        }
    }
}

#elif defined(HAVE_CONIO)
/* DOS/Windows with conio.h
 * Use putch() instead of putchar() for unbuffered output
 */

static void conio_puts(const char *s)
{
    while (*s) {
        putch(*s++);
    }
}

static void conio_backspace(void)
{
    putch('\b');
    putch(' ');
    putch('\b');
}

int read_line(char *buf, int max_len)
{
    int pos = 0;
    int len = 0;
    int c;
    
    history_pos = history_count;
    
    while (1) {
        c = getch();
        
        if (c == 0 || c == 0xE0) {
            int c2 = getch();
            if (c2 == 72) {  /* Up arrow */
                if (history_pos > 0) {
                    int i;
                    history_pos--;
                    for (i = 0; i < len; i++) conio_backspace();
                    strcpy(buf, history[history_pos]);
                    len = strlen(buf);
                    pos = len;
                    conio_puts(buf);
                }
            } else if (c2 == 80) {  /* Down arrow */
                if (history_pos < history_count - 1) {
                    int i;
                    history_pos++;
                    for (i = 0; i < len; i++) conio_backspace();
                    strcpy(buf, history[history_pos]);
                    len = strlen(buf);
                    pos = len;
                    conio_puts(buf);
                } else if (history_pos == history_count - 1) {
                    int i;
                    history_pos = history_count;
                    for (i = 0; i < len; i++) conio_backspace();
                    buf[0] = '\0';
                    len = 0;
                    pos = 0;
                }
            } else if (c2 == 75) {  /* Left arrow */
                if (pos > 0) {
                    pos--;
                    putch('\b');
                }
            } else if (c2 == 77) {  /* Right arrow */
                if (pos < len) {
                    putch(buf[pos]);
                    pos++;
                }
            }
            continue;
        }
        
        if (c == '\r' || c == '\n') {
            putch('\r');
            putch('\n');
            buf[len] = '\0';
            
            if (len > 0 && (history_count == 0 || strcmp(buf, history[history_count - 1]) != 0)) {
                if (history_count < MAX_HISTORY) {
                    strcpy(history[history_count++], buf);
                } else {
                    int i;
                    for (i = 0; i < MAX_HISTORY - 1; i++) {
                        strcpy(history[i], history[i + 1]);
                    }
                    strcpy(history[MAX_HISTORY - 1], buf);
                }
            }
            
            return len;
        }
        
        if (c == 8 || c == 127) {  /* Backspace */
            if (pos > 0) {
                int i;
                for (i = pos - 1; i < len - 1; i++) {
                    buf[i] = buf[i + 1];
                }
                len--;
                pos--;
                putch('\b');
                for (i = pos; i < len; i++) putch(buf[i]);
                putch(' ');
                for (i = len + 1; i > pos; i--) putch('\b');
            }
            continue;
        }
        
        if (c == 3) {  /* Ctrl-C */
            conio_puts("^C\r\n");
            buf[0] = '\0';
            return 0;
        }
        
        if (c >= 32 && c < 127 && len < max_len - 1) {
            int i;
            for (i = len; i > pos; i--) {
                buf[i] = buf[i - 1];
            }
            buf[pos] = c;
            len++;
            pos++;
            
            for (i = pos - 1; i < len; i++) putch(buf[i]);
            for (i = len; i > pos; i--) putch('\b');
        }
    }
}

#else
/* Fallback: simple fgets */
int read_line(char *buf, int max_len)
{
    if (fgets(buf, max_len, stdin) == NULL) {
        return -1;
    }
    
    {
        int len = strlen(buf);
        if (len > 0 && buf[len - 1] == '\n') {
            buf[--len] = '\0';
        }
        
        if (len > 0 && (history_count == 0 || strcmp(buf, history[history_count - 1]) != 0)) {
            if (history_count < MAX_HISTORY) {
                strcpy(history[history_count++], buf);
            }
        }
        
        return len;
    }
}
#endif
