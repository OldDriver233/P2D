import numpy as np
import matplotlib.pyplot as plt
import redis
import itertools

r = redis.Redis(host='localhost', port=6379, db=0)

print("Connected to redis")
data_len = r.llen('us')
data_list = r.lrange('us', 0, -1)

calc = np.fromiter(map(lambda x: float(x), data_list), dtype=np.float64)
x = np.arange(0, data_len, 1) / (data_len - 1)
plt.plot(x, calc)
plt.show()