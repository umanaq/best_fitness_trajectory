# 编译器和编译选项
CC = gcc
CFLAGS = -Wall -Wextra -O2 -Iinc
LDFLAGS = -lm

# 目录结构
SRC_DIR = src
INC_DIR = inc
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

# 目标可执行文件
TARGET = $(BIN_DIR)/pso_spline

# 源文件
SRCS = $(wildcard $(SRC_DIR)/*.c)
# 对象文件（放入obj子目录）
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
HEADERS = $(wildcard $(INC_DIR)/*.h)

# 创建必要目录（如果不存在）
MAKE_DIRS = $(OBJ_DIR) $(BIN_DIR)
$(shell mkdir -p $(MAKE_DIRS))

# 默认目标
all: $(TARGET)

# 链接生成可执行文件
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 编译规则
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# 清理规则
clean:
	rm -rf $(BUILD_DIR)

# 伪目标声明
.PHONY: all clean