#if defined(linux)

#include <limits.h>
#include <stdlib.h>

#include <string.h> /* strlen */

extern char *f2cstring(char* fstring, size_t len);
extern void padfstring(char *dest, const char *src, unsigned int len);

#define fullpathpp fullpathpp_
#define splitpathpp splitpathpp_
#define setenvpp setenvpp_
#define peekcharpp peekcharpp_


#if defined(__cplusplus)
extern "C" {
#endif

  int fullpathpp(char *name, char *path, unsigned int name_len, unsigned int path_len);

  int splitpathpp(char *path, char* drive, char* dir, char* name, char* ext,
		  unsigned int path_len, unsigned int drive_len, 
		  unsigned int dir_len, unsigned int name_len, 
		  unsigned int ext_len);
       
  int setenvpp(char *name, char *value, unsigned int name_len, unsigned int value_len);
  
  int peekcharpp(void);

#if defined(__cplusplus)
} 
#endif

       
int
fullpathpp(char *name, char *path, unsigned int name_len, unsigned int path_len)
{
  char *n;
  char *ptr;
  char buffer[PATH_MAX];

  n = f2cstring(name, (size_t) name_len);

  ptr = realpath(n, buffer);
  if (ptr == NULL) {
    // failure
    free(n);
    return 0;
  }
  else {
    padfstring(path, buffer, path_len);
    free(n);
    return strlen(buffer);
  }
}

int
splitpathpp(char *path, char* drive, char* dir, char* name, char* ext,
	    unsigned int path_len, unsigned int drive_len, 
	    unsigned int dir_len, unsigned int name_len, 
	    unsigned int ext_len)
{
  int i;
  size_t dot = 0;
  size_t slash = 0;
  int slash_found = 0;
  char *cpath = NULL;
  size_t plen = 0;

  /* linux has no drives */
  padfstring(drive, "", drive_len);

  if ((cpath = f2cstring(path, (size_t) path_len))) {
    plen = strlen(cpath);
  }

  if (plen == 0) {
    padfstring(ext, "", ext_len);      
    padfstring(name, "", name_len);      
    padfstring(dir, "", dir_len);
    free(cpath);
    return 0;
  }

  /* find last '.' and last '/' */
  for (i = plen - 1; i >= 0; --i) {
    if ((cpath[i] == '.') && dot == 0) {
      /* last '.' found */
      dot = i;
    }
    if (cpath[i] == '/') {
      /* last '/' found */
      slash_found = 1;
      slash = i;
      break;
    }
  }
  
  /* copy ext */
  if (dot) {
    /* ext found - '.' is included */
    padfstring(ext, &cpath[dot], ext_len);
    cpath[dot] = '\0';
  }
  else {
    /* no ext */
    padfstring(ext, "", ext_len);      
  }

  if (slash_found) {
    if (slash != plen - 1) {
      /* copy name */
      padfstring(name, &cpath[slash+1], name_len);
      cpath[slash+1] = '\0';
    }
    else {
      /* no name */
      padfstring(name, "", name_len);
    }

    /* copy dir  '/' is included */
    padfstring(dir, cpath, dir_len);    
  }
  else {
    /* copy name */
    padfstring(name, cpath, name_len);

    /* no dir */
    padfstring(dir, "", name_len);
  }
  
  plen = strlen(dir);
  free(cpath);
  return plen;
}

int
setenvpp(char *name, char *value, unsigned int name_len, unsigned int value_len)
{
  char *cname, *cvalue;
  int result = 0; /* failure */

  cname = f2cstring(name, (size_t) name_len);
  cvalue = f2cstring(value, (size_t) value_len);

  if (cname != NULL && cvalue != NULL) {
    if (setenv(cname, cvalue, 1) == 0) {
      result = 1; /* success */
    }
  }
  
  free(cname);
  free(cvalue);
  return result;
}

int
peekcharpp(void)
{
  return 0; /* ignore for now */
}


#endif /* defined(linux) */
