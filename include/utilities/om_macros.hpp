#if !defined(OM_MACROS)
#define OM_MACROS

#define OM_ASSERT(condition, message)                                          \
  do {                                                                         \
    if (!(condition)) {                                                        \
      std::cerr << "Assertion `" #condition "` failed in " << __FILE__         \
                << " line " << __LINE__ << ": " << message << std::endl;       \
      std::terminate();                                                        \
    }                                                                          \
  } while (false)

#endif /// OM_MACROS