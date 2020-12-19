#define DASH_LIST_LENGTH 2
#define DOTTED_LIST_LENGTH 2
#define DOT_DASHED_LIST_LENGTH 4
#define SHORT_DASHED_LIST_LENGTH 2
#define LONG_DASHED_LIST_LENGTH 2
#define ODD_DASHED_LIST_LENGTH 3

#define DASH_TOTAL_SIZE 20
#define DOTTED_TOTAL_SIZE 8
#define DOT_DASHED_TOTAL_SIZE 38
#define SHORT_DASHED_TOTAL_SIZE 20
#define LONG_DASHED_TOTAL_SIZE 34
#define ODD_DASHED_TOTAL_SIZE 24

	static int dash_total_size[] = {
		DASH_TOTAL_SIZE,
		DOTTED_TOTAL_SIZE,
		DOT_DASHED_TOTAL_SIZE,
		SHORT_DASHED_TOTAL_SIZE,
		LONG_DASHED_TOTAL_SIZE,
		ODD_DASHED_TOTAL_SIZE
	};

	static int dash_list_length[] = {
		DASH_LIST_LENGTH,
		DOTTED_LIST_LENGTH,
		DOT_DASHED_LIST_LENGTH,
		SHORT_DASHED_LIST_LENGTH,
		LONG_DASHED_LIST_LENGTH,
		ODD_DASHED_LIST_LENGTH
	};

	static char dash[DASH_LIST_LENGTH] = {10, 10};
	static char dotted[DOTTED_LIST_LENGTH] = {6, 2};
	static char dot_dashed[DOT_DASHED_LIST_LENGTH] = {26, 4, 4, 4};
	static char short_dashed[SHORT_DASHED_LIST_LENGTH] = {4, 16};
	static char long_dashed[LONG_DASHED_LIST_LENGTH] = {26, 8};
	static char odd_dashed[ODD_DASHED_LIST_LENGTH] = {4, 8, 12};
	static char *dash_list[] = {
		dash,
		dotted,
		dot_dashed,
		short_dashed,
		long_dashed,
		odd_dashed
	};
