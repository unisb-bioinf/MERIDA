/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef GENETRAIL2_CONFIG_MACROS_H
#define GENETRAIL2_CONFIG_MACROS_H

#define GT2_EXPORT          __attribute__((visibility ("default")))
#define GT2_LOCAL           __attribute__((visibility ("hidden")))
#define GT2_EXTERN_VARIABLE extern __attribute__((visibility ("default")))

#endif //GENETRAIL2_CONFIG_MACROS_H
