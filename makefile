# �������ͱ���ѡ��
CC = gcc
CFLAGS = -Wall -Wextra -O2 -Iinc
LDFLAGS = -lm

# Ŀ¼�ṹ
SRC_DIR = src
INC_DIR = inc
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

# Ŀ���ִ���ļ�
TARGET = $(BIN_DIR)/pso_spline

# Դ�ļ�
SRCS = $(wildcard $(SRC_DIR)/*.c)
# �����ļ�������obj��Ŀ¼��
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
HEADERS = $(wildcard $(INC_DIR)/*.h)

# ������ҪĿ¼����������ڣ�
MAKE_DIRS = $(OBJ_DIR) $(BIN_DIR)
$(shell mkdir -p $(MAKE_DIRS))

# Ĭ��Ŀ��
all: $(TARGET)

# �������ɿ�ִ���ļ�
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# �������
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# �������
clean:
	rm -rf $(BUILD_DIR)

# αĿ������
.PHONY: all clean